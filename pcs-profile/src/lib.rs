//! PCS-only profiling harness: the guest cost-model metrics for each multilinear PCS
//! (Basefold, Brakedown, and future schemes), measured in ISOLATION from any AIR.
//!
//! Two cost centers, matching the zkVM chip split (see `PCS_PROFILER_PLAN.md`):
//!   - **HASH**  — Poseidon2 permutations (the hash precompile), counted by role
//!                 (leaf / compress / sponge) via `p3-stark-stats`.
//!   - **ARITH** — ed25519 field ops (modmul / modadd / inverse chips): `mul / sq / inv / add`.
//!   - **SIZE**  — serialized opening-proof bytes (guest input).
//!
//! The hash runs over **BabyBear** (a different field than the ed25519 data), so the ed25519
//! counters measure ONLY the non-hash arithmetic — exactly the arithmetic-chip load. This mirrors
//! the guest `bb-hash` config. Build with `--features count-ops` (default) for the arithmetic numbers.

use p3_baby_bear::{BabyBear, Poseidon2BabyBear, default_babybear_poseidon2_16};
use p3_challenger::DuplexChallengerCF;
use p3_ed25519::Ed25519;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_stark_stats::{CountingPerm, HashCounts};
use p3_symmetric::{PaddingFreeSponge, SplitFieldHasher, TruncatedPermutation};

pub mod basefold;
pub mod brakedown;

// ---- BabyBear hash stack (role-counted), shared by every PCS ----

pub const BB_WIDTH: usize = 16;
pub const BB_RATE: usize = 8;
/// Merkle digest = 8 BabyBear ≈ 248-bit (~124-bit collision resistance).
pub const BB_DIGEST: usize = 8;

pub type PermBB = Poseidon2BabyBear<BB_WIDTH>;
/// The BabyBear permutation wrapped so every invocation is counted into a shared role counter.
pub type CPerm = CountingPerm<PermBB>;
pub type InnerHash = PaddingFreeSponge<CPerm, BB_WIDTH, BB_RATE, BB_DIGEST>;
/// Leaf hasher: ed25519 elements → BabyBear digest via base-|BabyBear| limbs.
pub type LeafHash = SplitFieldHasher<Ed25519, BabyBear, InnerHash>;
pub type Compress = TruncatedPermutation<CPerm, 2, BB_DIGEST, BB_WIDTH>;
pub type Mmcs = MerkleTreeMmcs<Ed25519, BabyBear, LeafHash, Compress, 2, BB_DIGEST>;
/// Transcript challenger: sponges over BabyBear, observes/samples ed25519.
pub type Challenger = DuplexChallengerCF<Ed25519, BabyBear, CPerm, BB_WIDTH, BB_RATE>;

/// Build the BabyBear Merkle MMCS with leaf/compress permutations counted into `counts`.
pub fn build_mmcs(counts: &HashCounts) -> Mmcs {
    let leaf = LeafHash::new(InnerHash::new(CountingPerm::with_counter(
        default_babybear_poseidon2_16(),
        counts.leaf_hash.clone(),
    )));
    let compress = Compress::new(CountingPerm::with_counter(
        default_babybear_poseidon2_16(),
        counts.compress.clone(),
    ));
    Mmcs::new(leaf, compress, 0)
}

/// Build the cross-field duplex challenger with its sponge permutation counted into `counts`.
pub fn build_challenger(counts: &HashCounts) -> Challenger {
    Challenger::new(CountingPerm::with_counter(
        default_babybear_poseidon2_16(),
        counts.sponge.clone(),
    ))
}

// ---- ed25519 arithmetic counters (verify-side) ----

/// Reset the ed25519 op counters (call just before the timed/counted verify).
pub fn reset_arith() {
    #[cfg(feature = "count-ops")]
    p3_ed25519::counters::reset();
}

/// `(mul, sq, inv, add)` since the last [`reset_arith`]. Zeros without `count-ops`.
pub fn snapshot_arith() -> (u64, u64, u64, u64) {
    #[cfg(feature = "count-ops")]
    {
        let c = p3_ed25519::counters::snapshot();
        (c.mul, c.sq, c.inv, c.add)
    }
    #[cfg(not(feature = "count-ops"))]
    {
        (0, 0, 0, 0)
    }
}

// ---- metrics ----

/// One verify-side measurement: the cost of a single commit/open/verify at a random point.
#[derive(Clone, Debug)]
pub struct PcsMetrics {
    pub pcs: &'static str,
    /// Proximity regime the query count derives from (e.g. "δ/3", "δ/2", "list").
    pub regime: &'static str,
    pub num_vars: usize,
    pub q: usize,
    /// Actual soundness achieved = `q · −log₂(1−e)` bits (< target when `q` clamps to codeword).
    pub sec_bits: f64,
    pub proof_bytes: usize,
    pub hash_leaf: u64,
    pub hash_compress: u64,
    pub hash_sponge: u64,
    pub mul: u64,
    pub sq: u64,
    pub inv: u64,
    pub add: u64,
    /// Commit (encode + Merkle) wall-clock, µs. Meaningful only with `count-ops` OFF.
    pub commit_us: u128,
    /// Open (sumcheck/fold + query extraction) wall-clock, µs. `count-ops` OFF.
    pub open_us: u128,
    /// Verify wall-clock, µs. `count-ops` OFF.
    pub verify_us: u128,
}

impl PcsMetrics {
    #[must_use]
    pub fn hash_total(&self) -> u64 {
        self.hash_leaf + self.hash_compress + self.hash_sponge
    }

    #[must_use]
    pub fn header() -> String {
        format!(
            "{:<11} {:>5} {:>3} {:>6} {:>6} {:>11} {:>10} {:>10} {:>9} {:>8} {:>8} {:>10} {:>10} {:>6} {:>11}",
            "PCS", "reg", "nv", "q", "sec", "proof_B", "commit_us", "open_us", "verify_us", "h_tot",
            "h_leaf", "mul", "sq", "inv", "add"
        )
    }

    #[must_use]
    pub fn row(&self) -> String {
        format!(
            "{:<11} {:>5} {:>3} {:>6} {:>6.1} {:>11} {:>10} {:>10} {:>9} {:>8} {:>8} {:>10} {:>10} {:>6} {:>11}",
            self.pcs,
            self.regime,
            self.num_vars,
            self.q,
            self.sec_bits,
            self.proof_bytes,
            self.commit_us,
            self.open_us,
            self.verify_us,
            self.hash_total(),
            self.hash_leaf,
            self.mul,
            self.sq,
            self.inv,
            self.add,
        )
    }
}

/// A multilinear PCS that can be measured PCS-only. Implementors drive their *core* type
/// (commit → open at a random point → verify) and return the verify-side [`PcsMetrics`].
pub trait ProfileablePcs {
    fn name(&self) -> &'static str;
    /// Commit a random `2^num_vars`-multilinear, open at a random point, then VERIFY under the
    /// counters. `seed` fixes the random poly/point; `q` is the query count (clamped internally).
    fn run(&self, num_vars: usize, q: usize, seed: u64) -> PcsMetrics;
}
