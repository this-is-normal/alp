//! Cross-field Merkle leaf hasher: hash BIG-field (`SF`) data with a cheap
//! SMALL 32-bit field (`TF`) permutation.
//!
//! Each `SF` element (e.g. a 255-bit ed25519 value) is decomposed into base-`|TF|`
//! limbs (injectively — enough limbs that distinct `SF` elements never collide
//! before hashing), and the resulting `TF` limb stream is fed to an inner
//! `TF`-native hasher `H: CryptographicHasher<TF, [TF; OUT]>` (e.g. a Poseidon2
//! sponge over BabyBear).
//!
//! Motivation: verifying a Merkle proof over `SF` data inside a 32-bit zkVM is
//! dominated by the hash. A native `SF` (255-bit) Poseidon2 permutation is
//! ~hundreds of emulated `u64`-limb multiplications; a BabyBear permutation is
//! native `u32` arithmetic (and could later be a VM precompile). Hashing the
//! limb-decomposition with BabyBear-Poseidon2 keeps the committed data in `SF`
//! (opened leaves stay `SF`, so no reconstruction is needed on the arithmetic
//! side) while making every permutation ~100x cheaper.

use core::marker::PhantomData;

use p3_field::{
    PrimeField, PrimeField32, inject_field_order_num_limbs, injective_pack_bits,
    split_pf_to_packed_limbs,
};

use crate::CryptographicHasher;

/// A [`CryptographicHasher`] mapping big-field `SF` items into a small-field `TF`
/// digest by base-`|TF|` limb decomposition + an inner `TF`-native hasher `H`.
#[derive(Clone, Debug)]
pub struct SplitFieldHasher<SF, TF, H> {
    inner: H,
    /// Injective limb count: smallest `n` with `|TF|^n >= |SF|`.
    num_limbs: usize,
    _phantom: PhantomData<(SF, TF)>,
}

impl<SF: PrimeField, TF: PrimeField32, H> SplitFieldHasher<SF, TF, H> {
    /// Wrap an inner `TF`-native hasher. The injective limb count is derived from
    /// the two field orders so `hash_iter` is collision-free before the inner hash.
    pub fn new(inner: H) -> Self {
        Self {
            inner,
            num_limbs: inject_field_order_num_limbs::<SF, TF>(),
            _phantom: PhantomData,
        }
    }
}

impl<SF, TF, H, const OUT: usize> CryptographicHasher<SF, [TF; OUT]> for SplitFieldHasher<SF, TF, H>
where
    SF: PrimeField,
    TF: PrimeField32,
    H: CryptographicHasher<TF, [TF; OUT]>,
{
    fn hash_iter<I>(&self, input: I) -> [TF; OUT]
    where
        I: IntoIterator<Item = SF>,
    {
        let n = self.num_limbs;
        // Radix-2^b bit-shift decomposition (b = injective_pack_bits::<TF>()) — avoids the
        // BigUint base-|TF| DIVISION that dominates the in-zkVM split cost. Still injective:
        // n limbs * b bits >= bits(SF), each limb < 2^b < |TF|.
        let b = injective_pack_bits::<TF>();
        self.inner.hash_iter(
            input
                .into_iter()
                .flat_map(|x| split_pf_to_packed_limbs::<SF, TF>(x, n, b)),
        )
    }
}
