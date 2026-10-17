#![doc = include_str!("../README.md")]
#![no_std]

extern crate alloc;

mod hiding_mmcs;
mod merkle_tree;
mod mmcs;
mod pruning;

pub use hiding_mmcs::*;
pub use merkle_tree::MerkleTree;
pub use mmcs::{MerkleTreeError, MerkleTreeMmcs, PrunedBatchOpening, PrunedProofError};
pub use p3_symmetric::MerkleCap;
pub use pruning::PrunedMerklePaths;

/// Counts pruned-batch verifications and the query-paths they cover.
///
/// The recursion verifier reaches Merkle work through `verify_batch_pruned`,
/// not the per-path `verify_batch` — an amortised frontier walk over all
/// queries at once. Multiplying these counts by the measured per-query cost
/// (`alp-openvm/guests/microbench` mode 12/13) is what turns "Merkle is
/// expensive" into a share of the guest.
#[cfg(feature = "count-ops")]
pub mod pruned_counter {
    use core::sync::atomic::{AtomicU64, Ordering};
    pub static CALLS: AtomicU64 = AtomicU64::new(0);
    pub static QUERY_PATHS: AtomicU64 = AtomicU64::new(0);
    #[inline]
    pub fn on_verify(num_queries: usize) {
        CALLS.fetch_add(1, Ordering::Relaxed);
        QUERY_PATHS.fetch_add(num_queries as u64, Ordering::Relaxed);
    }
    pub fn reset() {
        CALLS.store(0, Ordering::Relaxed);
        QUERY_PATHS.store(0, Ordering::Relaxed);
    }
    #[must_use]
    pub fn snapshot() -> (u64, u64) {
        (CALLS.load(Ordering::Relaxed), QUERY_PATHS.load(Ordering::Relaxed))
    }
}
