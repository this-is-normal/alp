//! **BaseFold** over an ed25519-native random foldable code: no FFT, and a
//! polylogarithmic verifier.

use p3_basefold::{BasefoldMlPcs, BasefoldPcs, RandomFoldableCode};
use p3_util::log2_ceil_usize;

use crate::hash::{F, L, MyMmcs, build_mmcs};

/// Human-readable backend name, for measurement output.
pub const NAME: &str = "BaseFold";

pub type CodePcs = BasefoldMlPcs<F, MyMmcs, L>;

/// Foldable-code inverse rate.
///
/// **1/16, not 1/4.** A lower rate buys a larger list-decoding radius, and the
/// query count falls with it: 241 queries at 1/4 against 136 at 1/16, for the
/// same 100-bit target net of grinding. The prover encodes a 4× longer
/// codeword, which is cheap (encoding is additions plus one twist multiply per
/// butterfly), and the *verifier* — which is what a recursion circuit runs —
/// does 1.8× less work, since its cost is queries × Merkle-path hashing.
const INV_RATE: usize = 16;
/// ChaCha seed (portable across host and zkVM).
const CODE_SEED: u64 = 0xBA5E_F01D;
/// List-decoding proximity radius at [`INV_RATE`]: 0.215 at rate 1/4, 0.350 at 1/16.
const RADIUS: f64 = 0.350;
/// Proof-of-work bits before the query phase, matched to the Quasar backend so
/// the comparison stays like-for-like. Each bit buys a bit of query soundness,
/// so the count is derived from the reduced target `100 − POW_BITS`.
const POW_BITS: usize = 16;

/// Build one BaseFold PCS sized for a `2^log_height`-row, `width`-column trace.
pub fn build_pcs(log_height: usize, width: usize) -> CodePcs {
    // BaseFold commits the stacked polynomial directly.
    let num_variables = log_height + log2_ceil_usize(width);
    let code = RandomFoldableCode::<F>::new(num_variables, INV_RATE, CODE_SEED);
    // `q = ⌈(100 − POW_BITS) / −log₂(1 − δ)⌉`, clamped to the codeword so tiny
    // instances stay well-formed.
    let target = (100.0 - POW_BITS as f64).max(1.0);
    let queries = (target / -(1.0 - RADIUS).log2()).ceil() as usize;
    let queries = core::cmp::min(queries, INV_RATE << num_variables);
    CodePcs::new(BasefoldPcs::new_with_grinding(
        code,
        build_mmcs(),
        queries,
        POW_BITS,
    ))
}
