//! A Poseidon2 permutation that a zkVM guest can hand to a precompile with **no
//! representation conversion**.
//!
//! # The problem
//!
//! OpenVM's Poseidon2 chip takes canonical `u32` lanes; `BabyBear` stores
//! Montgomery form. A guest wrapper therefore converts all 16 lanes out and
//! back on every permutation — and in a Merkle path the next permutation
//! immediately converts them again. Measured in an rv32im guest: **549
//! instructions per wrapped permutation against 16 for the bare chip call**,
//! which came to ~41% of the cost of verifying a Merkle path.
//!
//! # The fix
//!
//! Let `M(x) = R·x` be multiplication by the Montgomery radix `R = 2^32 mod P`,
//! a field bijection. Define
//!
//! ```text
//!     P'(v) = M⁻¹( P( M(v) ) )
//! ```
//!
//! `M(v)` is exactly the field element whose *canonical* value is `v`'s
//! Montgomery word. So a guest computes `P'` by reading the raw words, calling
//! the chip, and writing the raw words back: **zero conversions**. A host
//! computes it as two elementwise constant multiplications around the ordinary
//! permutation, which stays vectorised.
//!
//! `P'` is `P` conjugated by a bijection, so it is a permutation of the state
//! space with the same collision and preimage resistance. It is a *different*
//! function from `P`, so both sides must use this type — a mismatch shows up
//! immediately as a failed verification, not as a subtle weakness.

use p3_field::{Field, PrimeCharacteristicRing, PrimeField32};
use p3_symmetric::{CryptographicPermutation, Permutation};

use crate::BabyBear;

/// `R = 2^32 mod P` as a field element.
#[inline]
#[must_use]
pub fn monty_radix() -> BabyBear {
    BabyBear::from_u32(((1u64 << 32) % u64::from(BabyBear::ORDER_U32)) as u32)
}

/// Wraps a permutation `P` as `P'(v) = R⁻¹·P(R·v)` — see the module docs.
///
/// Use [`Self::permute_raw`] on a zkVM guest, where the whole point is that the
/// raw Montgomery words *are* the permutation's inputs.
#[derive(Clone, Debug)]
pub struct MontyBypass<P> {
    inner: P,
    radix: BabyBear,
    radix_inv: BabyBear,
}

impl<P> MontyBypass<P> {
    /// # Panics
    ///
    /// Never: `R` is a unit, being a power of two in an odd-order field.
    pub fn new(inner: P) -> Self {
        let radix = monty_radix();
        Self {
            inner,
            radix,
            radix_inv: radix.try_inverse().expect("R is a unit"),
        }
    }
}

impl<P: Permutation<[BabyBear; WIDTH]>, const WIDTH: usize> Permutation<[BabyBear; WIDTH]>
    for MontyBypass<P>
{
    fn permute_mut(&self, state: &mut [BabyBear; WIDTH]) {
        for x in state.iter_mut() {
            *x *= self.radix;
        }
        self.inner.permute_mut(state);
        for x in state.iter_mut() {
            *x *= self.radix_inv;
        }
    }
}

impl<P: Permutation<[BabyBear; WIDTH]>, const WIDTH: usize>
    CryptographicPermutation<[BabyBear; WIDTH]> for MontyBypass<P>
{
}

/// The identity a guest relies on: `P'`'s input lanes, in canonical form, are
/// exactly the input's raw Montgomery words, and likewise for the output.
///
/// A guest does not call this — it reads the words and calls the chip — but the
/// equality is what makes the two sides agree, so it is stated (and tested)
/// here rather than left implicit in guest code.
#[must_use]
pub fn raw_words<const WIDTH: usize>(state: &[BabyBear; WIDTH]) -> [u32; WIDTH] {
    core::array::from_fn(|i| state[i].to_monty_raw())
}

/// Inverse of [`raw_words`].
#[must_use]
pub fn from_raw_words<const WIDTH: usize>(words: &[u32; WIDTH]) -> [BabyBear; WIDTH] {
    core::array::from_fn(|i| BabyBear::from_monty_raw(words[i]))
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128Plus;

    use super::*;
    use crate::{Poseidon2BabyBear, default_babybear_poseidon2_16};

    /// What a guest computes — reinterpret the Montgomery words as canonical
    /// inputs, permute, reinterpret back — must equal `MontyBypass`.
    #[test]
    fn bypass_matches_the_raw_word_reinterpretation() {
        let inner = default_babybear_poseidon2_16();
        let bypass = MontyBypass::new(inner.clone());
        let mut rng = Xoroshiro128Plus::seed_from_u64(9);
        for _ in 0..8 {
            let state: [BabyBear; 16] = core::array::from_fn(|_| {
                BabyBear::from_u32(rand::RngExt::random::<u32>(&mut rng) % BabyBear::ORDER_U32)
            });

            // Host formulation.
            let mut host = state;
            bypass.permute_mut(&mut host);

            // Guest formulation: the raw words are the permutation's lanes.
            let words = raw_words(&state);
            let as_values: [BabyBear; 16] = core::array::from_fn(|i| BabyBear::from_u32(words[i]));
            let permuted = inner.permute(as_values);
            let out_words: [u32; 16] = core::array::from_fn(|i| permuted[i].as_canonical_u32());
            let guest = from_raw_words(&out_words);

            assert_eq!(host, guest, "host and guest formulations must agree");
        }
    }

    /// A permutation, not a collapse: distinct inputs stay distinct.
    #[test]
    fn bypass_is_injective_on_samples() {
        let bypass = MontyBypass::new(default_babybear_poseidon2_16());
        let a: [BabyBear; 16] = core::array::from_fn(|i| BabyBear::from_u32(i as u32));
        let mut b = a;
        b[3] += BabyBear::ONE;
        let (mut pa, mut pb) = (a, b);
        bypass.permute_mut(&mut pa);
        bypass.permute_mut(&mut pb);
        assert_ne!(pa, pb);
    }

    /// Sanity: the wrapper really is a different function from the inner one.
    #[test]
    fn bypass_differs_from_the_plain_permutation() {
        let inner: Poseidon2BabyBear<16> = default_babybear_poseidon2_16();
        let bypass = MontyBypass::new(inner.clone());
        let state: [BabyBear; 16] = core::array::from_fn(|i| BabyBear::from_u32(i as u32 + 1));
        let mut plain = state;
        inner.permute_mut(&mut plain);
        let mut wrapped = state;
        bypass.permute_mut(&mut wrapped);
        assert_ne!(plain, wrapped);
    }
}
