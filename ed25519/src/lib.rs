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

/// Throwaway global op-counters for profiling (M4). Counts ed25519 field mul / add+sub.
#[cfg(feature = "count-ops")]
pub mod counters {
    use core::sync::atomic::{AtomicU64, Ordering};
    pub static MUL: AtomicU64 = AtomicU64::new(0);
    pub static ADD: AtomicU64 = AtomicU64::new(0);
    pub static PERM: AtomicU64 = AtomicU64::new(0);
    #[inline]
    pub fn on_mul() {
        MUL.fetch_add(1, Ordering::Relaxed);
    }
    #[inline]
    pub fn on_add() {
        ADD.fetch_add(1, Ordering::Relaxed);
    }
    /// Count one Poseidon2 permutation (a hash invocation). Locally the hash is over Ed25519 so its
    /// muls land in MUL too; this lets us subtract them to isolate the PCS arithmetic muls.
    #[inline]
    pub fn on_perm() {
        PERM.fetch_add(1, Ordering::Relaxed);
    }
    pub fn reset() {
        MUL.store(0, Ordering::Relaxed);
        ADD.store(0, Ordering::Relaxed);
        PERM.store(0, Ordering::Relaxed);
    }
    /// (mul, add+sub, poseidon2-perms) since the last reset.
    pub fn snapshot() -> (u64, u64, u64) {
        (MUL.load(Ordering::Relaxed), ADD.load(Ordering::Relaxed), PERM.load(Ordering::Relaxed))
    }
}
pub use poseidon2::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, ED25519_S_BOX_DEGREE,
    Poseidon2Ed25519, Poseidon2ExternalLayerEd25519, Poseidon2InternalLayerEd25519,
};
