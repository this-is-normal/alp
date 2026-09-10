//! A parallel proof-of-work search.
//!
//! `GrindingChallenger::grind`'s default implementation tests candidates in a
//! serial loop. That made proof of work **the single largest item in the
//! prover**: a Quasar proof runs six independent grinds (a column phase plus two
//! BaseFold opening phases, for each of the two committed traces), and at 16
//! bits each expects `2^16` candidates. Measured at 64 signatures, the six
//! grinds were ~2.6 s of a 4.71 s prove — 55% — and none of it used more than
//! one core, which capped the achievable speedup at about 1.8x however many
//! cores were available.
//!
//! Lives here rather than in `p3-challenger` so the core transcript crate keeps
//! its dependency set; `p3-quasar` reaches it through `p3-basefold`.

use p3_challenger::GrindingChallenger;
use p3_field::PrimeCharacteristicRing;
use p3_maybe_rayon::prelude::*;

/// Candidates tested per parallel window. Small enough that the final window
/// overshoots the found witness by only a few percent of the expected work,
/// large enough to keep every core busy.
const WINDOW: u64 = 4096;

/// The **smallest** witness `w` with `check_witness(bits, w)`, searched in
/// parallel.
///
/// Returning the smallest — not merely the first found — is what keeps the
/// search deterministic: the same transcript state always yields the same
/// witness, so proofs stay reproducible across thread counts and backends.
///
/// Leaves `challenger` in exactly the state the serial `grind` would: the
/// witness observed and its bits sampled. That is required, not cosmetic — the
/// verifier's `check_witness` performs the same absorption, and any mismatch
/// desyncs every subsequent challenge.
pub fn grind_parallel<Ch>(challenger: &mut Ch, bits: usize) -> Ch::Witness
where
    Ch: GrindingChallenger + Clone + Send + Sync,
    Ch::Witness: PrimeCharacteristicRing,
{
    if bits == 0 {
        return Ch::Witness::ZERO;
    }

    let mut base = 0u64;
    let witness = loop {
        // `min` rather than a short-circuiting search: it is the portable
        // choice across `p3-maybe-rayon`'s two backends, and it is what makes
        // the result the smallest candidate rather than whichever thread
        // finished first.
        let found = (base..base + WINDOW)
            .into_par_iter()
            .filter(|&i| {
                let mut probe = challenger.clone();
                probe.check_witness(bits, Ch::Witness::from_u64(i))
            })
            .min();
        if let Some(i) = found {
            break Ch::Witness::from_u64(i);
        }
        base += WINDOW;
    };

    // Advance the real transcript exactly as the verifier will.
    assert!(
        challenger.check_witness(bits, witness),
        "the witness found in parallel must check against the real challenger"
    );
    witness
}
