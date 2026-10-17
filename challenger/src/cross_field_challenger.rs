//! A duplex-sponge challenger whose transcript sponge lives in a SMALL 32-bit
//! field `PF` (e.g. BabyBear) but whose observe/sample interface is a BIG prime
//! field `F` (e.g. the 255-bit ed25519 base field).
//!
//! This is the mirror image of [`MultiField32Challenger`](crate::MultiField32Challenger)
//! (which sponges over a big field and samples a small one). Here:
//!
//! - **observe(`F`)**: decompose the big element into base-`|PF|` limbs
//!   (injective) and absorb them into the `PF` sponge — cheap (a split is bit
//!   shuffling).
//! - **sample() → `F`**: squeeze `num_limbs` `PF` cells and recombine them
//!   (Horner in base `|PF|`) into a big-field challenge. This reconstruction is
//!   the *only* place a big element is rebuilt from limbs — it runs
//!   `O(#challenges)` times (tens), never per Merkle leaf.
//! - **observe(`MerkleCap<F,[PF;N]>`) / `Hash<F,PF,N>`**: absorb the `PF` digest
//!   words natively (the digests are already small-field).
//!
//! Rationale: pairs with [`SplitFieldHasher`](p3_symmetric::SplitFieldHasher) so a
//! Plonky3 proof over `F` uses a cheap `PF`-Poseidon2 for all Merkle/Fiat-Shamir
//! hashing, which is ~100x cheaper per permutation than a native `F` permutation
//! inside a 32-bit zkVM verifier.

use alloc::vec::Vec;
use core::marker::PhantomData;

use p3_field::{
    BasedVectorSpace, PrimeField, PrimeField32, inject_field_order_num_limbs, injective_pack_bits,
    split_pf_to_packed_limbs,
};
use p3_symmetric::{CryptographicPermutation, Hash, MerkleCap};

use crate::{
    CanObserve, CanSample, CanSampleBits, CanSampleUniformBits, DuplexChallenger, FieldChallenger,
    GrindingChallenger, ResamplingError,
};

/// A [`DuplexChallenger`] over the small field `PF`, adapted to observe/sample the
/// big field `F` via base-`|PF|` limb (de)composition.
#[derive(Clone, Debug)]
pub struct DuplexChallengerCF<F, PF, P, const WIDTH: usize, const RATE: usize>
where
    PF: Clone,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    inner: DuplexChallenger<PF, P, WIDTH, RATE>,
    /// Base-`|PF|` limbs per big-field element (both split on observe and combine on sample).
    num_limbs: usize,
    _phantom: PhantomData<F>,
}

impl<F, PF, P, const WIDTH: usize, const RATE: usize> DuplexChallengerCF<F, PF, P, WIDTH, RATE>
where
    F: PrimeField,
    PF: PrimeField32,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    /// Build a fresh challenger over `perm` (a `PF` permutation). The limb count is
    /// derived from `|F|` and `|PF|` so observe is injective and sample near-uniform.
    pub fn new(perm: P) -> Self {
        Self {
            inner: DuplexChallenger::new(perm),
            num_limbs: inject_field_order_num_limbs::<F, PF>(),
            _phantom: PhantomData,
        }
    }

    /// Squeeze `num_limbs` `PF` cells and recombine (base-`|PF|` Horner) into one `F`.
    fn sample_big(&mut self) -> F {
        let base = F::from_u32(PF::ORDER_U32);
        let limbs: Vec<PF> = (0..self.num_limbs).map(|_| self.inner.sample()).collect();
        // limbs are little-endian (limb[0] least significant); fold from the top.
        let mut acc = F::ZERO;
        for limb in limbs.iter().rev() {
            acc = acc * base + F::from_u32(limb.as_canonical_u32());
        }
        acc
    }
}

impl<F, PF, P, const WIDTH: usize, const RATE: usize> CanObserve<F>
    for DuplexChallengerCF<F, PF, P, WIDTH, RATE>
where
    F: PrimeField,
    PF: PrimeField32,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    fn observe(&mut self, value: F) {
        // Radix-2^b bit-shift split (cheap in-zkVM; avoids BigUint division). Injective.
        let b = injective_pack_bits::<PF>();
        for limb in split_pf_to_packed_limbs::<F, PF>(value, self.num_limbs, b) {
            self.inner.observe(limb);
        }
    }
}

impl<F, PF, const N: usize, P, const WIDTH: usize, const RATE: usize> CanObserve<Hash<F, PF, N>>
    for DuplexChallengerCF<F, PF, P, WIDTH, RATE>
where
    F: PrimeField,
    PF: PrimeField32,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    fn observe(&mut self, values: Hash<F, PF, N>) {
        for w in values {
            self.inner.observe(w);
        }
    }
}

impl<F, PF, const N: usize, P, const WIDTH: usize, const RATE: usize>
    CanObserve<&MerkleCap<F, [PF; N]>> for DuplexChallengerCF<F, PF, P, WIDTH, RATE>
where
    F: PrimeField,
    PF: PrimeField32,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    fn observe(&mut self, cap: &MerkleCap<F, [PF; N]>) {
        for digest in cap.roots() {
            for &w in digest {
                self.inner.observe(w);
            }
        }
    }
}

impl<F, PF, const N: usize, P, const WIDTH: usize, const RATE: usize>
    CanObserve<MerkleCap<F, [PF; N]>> for DuplexChallengerCF<F, PF, P, WIDTH, RATE>
where
    F: PrimeField,
    PF: PrimeField32,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    fn observe(&mut self, cap: MerkleCap<F, [PF; N]>) {
        self.observe(&cap);
    }
}

impl<F, EF, PF, P, const WIDTH: usize, const RATE: usize> CanSample<EF>
    for DuplexChallengerCF<F, PF, P, WIDTH, RATE>
where
    F: PrimeField,
    EF: BasedVectorSpace<F>,
    PF: PrimeField32,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    fn sample(&mut self) -> EF {
        EF::from_basis_coefficients_fn(|_| self.sample_big())
    }
}

impl<F, PF, P, const WIDTH: usize, const RATE: usize> CanSampleBits<usize>
    for DuplexChallengerCF<F, PF, P, WIDTH, RATE>
where
    F: PrimeField,
    PF: PrimeField32,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    /// Query indices need only a handful of bits (`< log2(codeword) ≪ 31`), so one
    /// `PF` element supplies them; its low `bits` are statistically uniform.
    fn sample_bits(&mut self, bits: usize) -> usize {
        debug_assert!(bits < 31);
        let rand_pf: PF = self.inner.sample();
        (rand_pf.as_canonical_u32() as usize) & ((1 << bits) - 1)
    }
}

impl<F, PF, P, const WIDTH: usize, const RATE: usize> FieldChallenger<F>
    for DuplexChallengerCF<F, PF, P, WIDTH, RATE>
where
    F: PrimeField,
    PF: PrimeField32,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
}

impl<F, PF, P, const WIDTH: usize, const RATE: usize> CanSampleUniformBits<F>
    for DuplexChallengerCF<F, PF, P, WIDTH, RATE>
where
    F: PrimeField,
    PF: PrimeField32,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    fn sample_uniform_bits<const RESAMPLE: bool>(
        &mut self,
        bits: usize,
    ) -> Result<usize, ResamplingError> {
        Ok(self.sample_bits(bits))
    }
}

impl<F, PF, P, const WIDTH: usize, const RATE: usize> GrindingChallenger
    for DuplexChallengerCF<F, PF, P, WIDTH, RATE>
where
    F: PrimeField,
    PF: PrimeField32,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    type Witness = F;

    /// Brute-force PoW search (the Brakedown path always calls with `bits = 0`, so
    /// `check_witness` short-circuits and this returns immediately).
    fn grind(&mut self, bits: usize) -> Self::Witness {
        let mut i: u64 = 0;
        let witness = loop {
            let candidate = F::from_u64(i);
            if self.clone().check_witness(bits, candidate) {
                break candidate;
            }
            i += 1;
        };
        assert!(self.check_witness(bits, witness));
        witness
    }
}
