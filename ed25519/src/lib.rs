//! The ed25519 base field, defined as `F_p` where `p = 2^255 - 19`.
//!
//! This is a 255-bit prime field implemented in the same style as [`p3_bn254`], using a
//! 4-limb little-endian Montgomery representation with respect to `R = 2^256`.
//!
//! Unlike bn254, this field is **not** two-adic-friendly (`p - 1 = 2^2 * (2^253 - 5)` has
//! 2-adicity only 2), so it does **not** implement `TwoAdicField`. It is also not
//! `PrimeField64`/`PrimeField32`.
#![no_std]

extern crate alloc;

mod ed25519;
mod helpers;
mod poseidon2;

pub use ed25519::*;

#[cfg(feature = "count-ops")]
pub mod counters {
    //! Throwaway op counters for profiling verify. Separates ed25519 field ops
    //! (add/sub, mul, square, inverse) from hash permutations. In the BabyBear-hash
    //! config the hashing runs in BabyBear, so `HASH_COUNT` (BabyBear Poseidon2 perms)
    //! and the ed25519 field counters are cleanly disjoint.
    use core::sync::atomic::{AtomicU64, Ordering};
    pub static ADD_TOTAL: AtomicU64 = AtomicU64::new(0);  // Add + Sub
    pub static MUL_TOTAL: AtomicU64 = AtomicU64::new(0);  // Mul (excludes square())
    pub static SQ_TOTAL: AtomicU64 = AtomicU64::new(0);   // square()
    pub static INV_TOTAL: AtomicU64 = AtomicU64::new(0);  // try_inverse (each ~ exp by p-2)
    pub static PERM_COUNT: AtomicU64 = AtomicU64::new(0); // ed25519 Poseidon2 (unused in bb config)
    pub static HASH_COUNT: AtomicU64 = AtomicU64::new(0); // BabyBear Poseidon2 permutations
    #[inline] pub fn on_add() { ADD_TOTAL.fetch_add(1, Ordering::Relaxed); }
    #[inline] pub fn on_mul() { MUL_TOTAL.fetch_add(1, Ordering::Relaxed); }
    #[inline] pub fn on_sq() { SQ_TOTAL.fetch_add(1, Ordering::Relaxed); }
    #[inline] pub fn on_inv() { INV_TOTAL.fetch_add(1, Ordering::Relaxed); }
    #[inline] pub fn on_perm() { PERM_COUNT.fetch_add(1, Ordering::Relaxed); }
    #[inline] pub fn on_hash() { HASH_COUNT.fetch_add(1, Ordering::Relaxed); }
    pub fn reset() {
        for a in [&ADD_TOTAL, &MUL_TOTAL, &SQ_TOTAL, &INV_TOTAL, &PERM_COUNT, &HASH_COUNT] {
            a.store(0, Ordering::Relaxed);
        }
    }
    #[derive(Clone, Copy, Debug, Default)]
    pub struct Counts { pub add: u64, pub mul: u64, pub sq: u64, pub inv: u64, pub hash: u64 }
    pub fn snapshot() -> Counts {
        Counts {
            add: ADD_TOTAL.load(Ordering::Relaxed),
            mul: MUL_TOTAL.load(Ordering::Relaxed),
            sq: SQ_TOTAL.load(Ordering::Relaxed),
            inv: INV_TOTAL.load(Ordering::Relaxed),
            hash: HASH_COUNT.load(Ordering::Relaxed),
        }
    }
}
pub use poseidon2::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, ED25519_S_BOX_DEGREE,
    Poseidon2Ed25519, Poseidon2ExternalLayerEd25519, Poseidon2InternalLayerEd25519,
};
