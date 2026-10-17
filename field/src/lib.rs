#![doc = include_str!("../README.md")]
#![no_std]

extern crate alloc;

mod array;
mod batch_inverse;
pub mod coset;
mod dup;
pub mod exponentiation;
pub mod extension;
mod field;
mod helpers;
pub mod integers;
pub mod op_assign_macros;
mod packed;
mod sqrt;

pub use array::*;
pub use batch_inverse::*;
pub use dup::Dup;
pub use field::*;
pub use helpers::*;

/// Counts calls to [`split_pf_to_packed_limbs_into`], the cross-field limb
/// decomposition every hashed element goes through. Separate from the ed25519
/// op counters because it is neither field arithmetic nor a hash permutation —
/// it is the glue between them, and in a 32-bit zkVM guest it used to cost more
/// than either.
#[cfg(feature = "count-ops")]
pub mod split_counter {
    use core::sync::atomic::{AtomicU64, Ordering};
    pub static SPLIT_COUNT: AtomicU64 = AtomicU64::new(0);
    #[inline]
    pub fn on_split() {
        SPLIT_COUNT.fetch_add(1, Ordering::Relaxed);
    }
    pub fn reset() {
        SPLIT_COUNT.store(0, Ordering::Relaxed);
    }
    #[must_use]
    pub fn snapshot() -> u64 {
        SPLIT_COUNT.load(Ordering::Relaxed)
    }
}
pub use packed::*;
pub use sqrt::{tonelli_shanks, tonelli_shanks_two_adic};
