//! A duplex-sponge challenger that works over a general [`PrimeField`], including
//! fields that are **not** [`PrimeField64`](p3_field::PrimeField64) (e.g. the
//! 255-bit ed25519 base field).
//!
//! The core duplex sponge in [`DuplexChallenger`] already works for any
//! `F: Field` for the observe/sample operations. However, the challenger traits
//! that the Brakedown PCS needs on top of that — `CanSampleBits<usize>`,
//! `GrindingChallenger`, `CanSampleUniformBits`, and the `FieldChallenger`
//! marker — are all gated on `F: PrimeField64` on the bare [`DuplexChallenger`]
//! because they reduce a sampled field element to bits via `as_canonical_u64`.
//!
//! We cannot add blanket `PrimeField` impls of those traits directly on
//! [`DuplexChallenger`] without overlapping the existing `PrimeField64` impls.
//! [`DuplexChallengerPF`] is a thin newtype wrapper that provides them for any
//! `F: PrimeField`, reducing sampled elements to bits via the field's canonical
//! `BigUint` representation instead of `as_canonical_u64`.

use p3_field::PrimeField;
use p3_symmetric::CryptographicPermutation;

use crate::{
    CanObserve, CanSample, CanSampleBits, CanSampleUniformBits, DuplexChallenger, FieldChallenger,
    GrindingChallenger, ResamplingError,
};

/// A [`DuplexChallenger`] wrapper usable over any [`PrimeField`], not just
/// [`PrimeField64`](p3_field::PrimeField64).
///
/// It delegates observation and sampling to the inner [`DuplexChallenger`] (whose
/// sponge core is field-agnostic), and supplies bit-sampling / grinding /
/// uniform-bit-sampling / the [`FieldChallenger`] marker via the field's
/// canonical `BigUint` representation.
#[derive(Clone, Debug)]
pub struct DuplexChallengerPF<F, P, const WIDTH: usize, const RATE: usize>(
    pub DuplexChallenger<F, P, WIDTH, RATE>,
)
where
    F: Clone,
    P: CryptographicPermutation<[F; WIDTH]>;

impl<F, P, const WIDTH: usize, const RATE: usize> DuplexChallengerPF<F, P, WIDTH, RATE>
where
    F: Copy + Default,
    P: CryptographicPermutation<[F; WIDTH]>,
{
    /// Build a fresh challenger from the given permutation, zero-initialising the
    /// sponge state (identical to [`DuplexChallenger::new`]).
    pub fn new(perm: P) -> Self {
        Self(DuplexChallenger::new(perm))
    }
}

// --- Delegation: observe / sample forward to the inner duplex sponge. ---
//
// These single generic impls give us `CanObserve` for `F`, `[F; N]`,
// `Hash<F, F, N>`, and the Merkle commitment types, plus `CanSample` for `F` and
// any extension field — for free, since the inner `DuplexChallenger` already
// implements all of them for `F: Field`.

impl<F, P, const WIDTH: usize, const RATE: usize, T> CanObserve<T>
    for DuplexChallengerPF<F, P, WIDTH, RATE>
where
    F: Clone,
    P: CryptographicPermutation<[F; WIDTH]>,
    DuplexChallenger<F, P, WIDTH, RATE>: CanObserve<T>,
{
    #[inline]
    fn observe(&mut self, value: T) {
        self.0.observe(value);
    }

    #[inline]
    fn observe_slice(&mut self, values: &[T])
    where
        T: Clone,
    {
        self.0.observe_slice(values);
    }
}

impl<F, P, const WIDTH: usize, const RATE: usize, EF> CanSample<EF>
    for DuplexChallengerPF<F, P, WIDTH, RATE>
where
    F: Clone,
    P: CryptographicPermutation<[F; WIDTH]>,
    DuplexChallenger<F, P, WIDTH, RATE>: CanSample<EF>,
{
    #[inline]
    fn sample(&mut self) -> EF {
        self.0.sample()
    }
}

/// Return the least-significant 64-bit limb of a canonical `BigUint`.
///
/// For a field element `f`, `low_u64(f.as_canonical_biguint())` yields the low
/// 64 bits of `f`'s canonical integer representative. For a ~255-bit field this
/// is essentially uniform for the small bit widths (< 63) used by PCS queries.
#[inline]
fn low_u64(x: &num_bigint::BigUint) -> u64 {
    x.iter_u64_digits().next().unwrap_or(0)
}

impl<F, P, const WIDTH: usize, const RATE: usize> CanSampleBits<usize>
    for DuplexChallengerPF<F, P, WIDTH, RATE>
where
    F: PrimeField,
    P: CryptographicPermutation<[F; WIDTH]>,
{
    /// Sample a field element and take its low `bits` bits.
    ///
    /// For a large (~255-bit) prime field the low bits of a uniformly sampled
    /// element are statistically indistinguishable from uniform for the small
    /// `bits` values (< 63) used here.
    fn sample_bits(&mut self, bits: usize) -> usize {
        debug_assert!(bits < 63);
        let rand_f: F = self.0.sample();
        let low = low_u64(&rand_f.as_canonical_biguint());
        (low & ((1u64 << bits) - 1)) as usize
    }
}

impl<F, P, const WIDTH: usize, const RATE: usize> FieldChallenger<F>
    for DuplexChallengerPF<F, P, WIDTH, RATE>
where
    F: PrimeField,
    P: CryptographicPermutation<[F; WIDTH]>,
{
    // `FieldChallenger<F>` is a marker whose default methods rely on the
    // `CanObserve<F> + CanSample<F> + CanSampleBits<usize> + Sync` supertraits,
    // all of which are satisfied above.
}

impl<F, P, const WIDTH: usize, const RATE: usize> GrindingChallenger
    for DuplexChallengerPF<F, P, WIDTH, RATE>
where
    F: PrimeField,
    P: CryptographicPermutation<[F; WIDTH]>,
{
    type Witness = F;

    /// Brute-force search for a proof-of-work witness.
    ///
    /// Iterates candidate witnesses `0, 1, 2, …` (as field elements), returning
    /// the first that satisfies [`check_witness`](GrindingChallenger::check_witness).
    ///
    /// Note: the Brakedown path always calls this with `pow_bits = 0`, so
    /// `check_witness` short-circuits to `true` on the first candidate and no real
    /// work is performed. The loop is written to be correct for any `bits`, and —
    /// unlike the `PrimeField64` `DuplexChallenger` impl — it does not bound the
    /// search by `ORDER_U64` (which is undefined for a 255-bit field).
    fn grind(&mut self, bits: usize) -> Self::Witness {
        let mut i: u64 = 0;
        let witness = loop {
            let candidate = F::from_u64(i);
            // Test on a clone so the real transcript state is untouched until we
            // commit to the winning witness below.
            if self.clone().check_witness(bits, candidate) {
                break candidate;
            }
            i += 1;
        };
        // Re-run on `self` to advance the real transcript and re-confirm.
        assert!(self.check_witness(bits, witness));
        witness
    }
}

impl<F, P, const WIDTH: usize, const RATE: usize> CanSampleUniformBits<F>
    for DuplexChallengerPF<F, P, WIDTH, RATE>
where
    F: PrimeField,
    P: CryptographicPermutation<[F; WIDTH]>,
{
    /// For a ~255-bit field, the low `bits` of a sampled element are already
    /// (statistically) uniform for the small `bits` used here, so uniform bit
    /// sampling coincides with [`CanSampleBits::sample_bits`] and never needs to
    /// resample. Hence this is infallible regardless of `RESAMPLE`.
    fn sample_uniform_bits<const RESAMPLE: bool>(
        &mut self,
        bits: usize,
    ) -> Result<usize, ResamplingError> {
        Ok(self.sample_bits(bits))
    }
}

#[cfg(test)]
mod tests {
    use p3_baby_bear::BabyBear;
    use p3_field::PrimeCharacteristicRing;
    use p3_symmetric::{CryptographicPermutation, Permutation};

    use super::*;
    use crate::GrindingChallenger;

    const WIDTH: usize = 8;
    const RATE: usize = 4;

    type BB = BabyBear;

    /// A trivial (non-cryptographic) test permutation that reverses the state.
    /// Sufficient to exercise the wrapper's generic (non-`PrimeField64`) code path.
    #[derive(Clone)]
    struct RevPermutation;

    impl<F: Clone> Permutation<[F; WIDTH]> for RevPermutation {
        fn permute_mut(&self, input: &mut [F; WIDTH]) {
            input.reverse();
        }
    }
    impl<F: Clone> CryptographicPermutation<[F; WIDTH]> for RevPermutation {}

    type Chal = DuplexChallengerPF<BB, RevPermutation, WIDTH, RATE>;

    #[test]
    fn sample_bits_in_range() {
        let mut chal = Chal::new(RevPermutation);
        for i in 0..20u8 {
            chal.observe(BB::from_u8(i));
        }
        for bits in 0..10usize {
            for _ in 0..50 {
                let v = chal.sample_bits(bits);
                assert!(v < (1usize << bits), "sample_bits({bits}) = {v} out of range");
            }
        }
    }

    #[test]
    fn observe_sample_deterministic() {
        let mut a = Chal::new(RevPermutation);
        let mut b = Chal::new(RevPermutation);
        for i in 0..15u8 {
            a.observe(BB::from_u8(i));
            b.observe(BB::from_u8(i));
        }
        // Identical observation streams produce identical samples.
        for _ in 0..10 {
            let sa: BB = a.sample();
            let sb: BB = b.sample();
            assert_eq!(sa, sb);
        }
        // And identical bit samples.
        for bits in 1..8usize {
            assert_eq!(a.sample_bits(bits), b.sample_bits(bits));
        }
    }

    #[test]
    fn grind_zero_bits_is_trivial() {
        // With 0 pow bits (the Brakedown path — the only way `grind` is invoked in
        // practice), `grind` returns immediately and `check_witness` accepts.
        let mut chal = Chal::new(RevPermutation);
        let w = chal.grind(0);
        assert!(chal.check_witness(0, w));
    }
}
