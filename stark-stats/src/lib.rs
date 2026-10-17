//! Profiling / introspection helpers for the multilinear (`multi-stark`) +
//! Brakedown proving path.
//!
//! Two independent tools, both designed to stay **out of any timed path**
//! (counts are deterministic — measure them in a separate run; never instrument
//! the run you are timing):
//!
//! * **Counting wrappers** ([`CountingPerm`], [`CountingHasher`],
//!   [`CountingCompress`]) — newtypes over the `p3-symmetric` hash primitives
//!   that delegate to the inner value and bump a shared [`OpCounter`]. Plug them
//!   into a config's hash roles (leaf hasher / node compressor / challenger
//!   permutation) to count exactly how many of each the prover and verifier do.
//! * **Static reports** ([`CircuitShape`], [`PcsInfo`]) + [`serialized_size`] —
//!   cheap introspection of the circuit shape, the committed-polynomial /
//!   commitment parameters, and the byte size of a proof.
//!
//! Because the wrappers are generic over the `p3-symmetric` traits and hold only
//! an `Arc<AtomicUsize>`, a build that uses the raw (un-wrapped) primitives is
//! unaffected — the instrumentation exists only when you choose the wrapped
//! types.

#![no_std]

extern crate alloc;

use alloc::sync::Arc;
use alloc::vec::Vec;
use core::fmt;
use core::sync::atomic::{AtomicUsize, Ordering};

use p3_symmetric::{CryptographicHasher, CryptographicPermutation, Permutation, PseudoCompressionFunction};
use serde::Serialize;

// ---------------------------------------------------------------------------
// Counter
// ---------------------------------------------------------------------------

/// A shared, resettable operation counter.
///
/// Cloning shares the same underlying atomic, so all clones of a wrapped
/// primitive (the prover/verifier clone hashers freely) accumulate into one
/// count. Uses `Relaxed` ordering — we only need the final total, and the point
/// is to count in a *separate, untimed* run, so contention cost is irrelevant.
#[derive(Clone, Default)]
pub struct OpCounter(Arc<AtomicUsize>);

impl OpCounter {
    /// A fresh counter at zero.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// The current count.
    #[must_use]
    pub fn count(&self) -> usize {
        self.0.load(Ordering::Relaxed)
    }

    /// Reset to zero (e.g. between the prove and verify phases).
    pub fn reset(&self) {
        self.0.store(0, Ordering::Relaxed);
    }

    #[inline]
    fn inc(&self) {
        self.0.fetch_add(1, Ordering::Relaxed);
    }
}

// ---------------------------------------------------------------------------
// Counting wrappers
// ---------------------------------------------------------------------------

/// Wraps a [`CryptographicPermutation`], counting each `permute`/`permute_mut`.
#[derive(Clone)]
pub struct CountingPerm<P> {
    inner: P,
    counter: OpCounter,
}

impl<P> CountingPerm<P> {
    /// Wrap `inner` with a fresh counter.
    pub fn new(inner: P) -> Self {
        Self { inner, counter: OpCounter::new() }
    }

    /// Wrap `inner`, sharing an existing counter handle (keep a clone to read it).
    pub const fn with_counter(inner: P, counter: OpCounter) -> Self {
        Self { inner, counter }
    }

    /// A clone of the counter handle.
    pub fn counter(&self) -> OpCounter {
        self.counter.clone()
    }
}

impl<T: Clone, P: Permutation<T>> Permutation<T> for CountingPerm<P> {
    fn permute(&self, input: T) -> T {
        self.counter.inc();
        self.inner.permute(input)
    }

    fn permute_mut(&self, input: &mut T) {
        self.counter.inc();
        self.inner.permute_mut(input);
    }
}

impl<T: Clone, P: CryptographicPermutation<T>> CryptographicPermutation<T> for CountingPerm<P> {}

/// Wraps a [`CryptographicHasher`], counting each hash invocation.
#[derive(Clone)]
pub struct CountingHasher<H> {
    inner: H,
    counter: OpCounter,
}

impl<H> CountingHasher<H> {
    pub fn new(inner: H) -> Self {
        Self { inner, counter: OpCounter::new() }
    }

    pub const fn with_counter(inner: H, counter: OpCounter) -> Self {
        Self { inner, counter }
    }

    pub fn counter(&self) -> OpCounter {
        self.counter.clone()
    }
}

impl<Item: Clone, Out, H: CryptographicHasher<Item, Out>> CryptographicHasher<Item, Out>
    for CountingHasher<H>
{
    fn hash_iter<I>(&self, input: I) -> Out
    where
        I: IntoIterator<Item = Item>,
    {
        // Count once per logical hash; the trait's other methods default to this.
        self.counter.inc();
        self.inner.hash_iter(input)
    }
}

/// Wraps a [`PseudoCompressionFunction`], counting each 2-to-1 compression.
#[derive(Clone)]
pub struct CountingCompress<C> {
    inner: C,
    counter: OpCounter,
}

impl<C> CountingCompress<C> {
    pub fn new(inner: C) -> Self {
        Self { inner, counter: OpCounter::new() }
    }

    pub const fn with_counter(inner: C, counter: OpCounter) -> Self {
        Self { inner, counter }
    }

    pub fn counter(&self) -> OpCounter {
        self.counter.clone()
    }
}

impl<T, const N: usize, C: PseudoCompressionFunction<T, N>> PseudoCompressionFunction<T, N>
    for CountingCompress<C>
{
    fn compress(&self, input: [T; N]) -> T {
        self.counter.inc();
        self.inner.compress(input)
    }
}

/// A per-role hash-op tally: leaf hashing, node compression, and the challenger
/// sponge permutation. Share these three counters into the wrapped roles of an
/// instrumented config, then read them after prove / verify.
#[derive(Clone, Default)]
pub struct HashCounts {
    pub leaf_hash: OpCounter,
    pub compress: OpCounter,
    pub sponge: OpCounter,
}

impl HashCounts {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset all three (e.g. between prove and verify).
    pub fn reset(&self) {
        self.leaf_hash.reset();
        self.compress.reset();
        self.sponge.reset();
    }

    /// `(leaf_hash, compress, sponge)` current counts.
    #[must_use]
    pub fn snapshot(&self) -> (usize, usize, usize) {
        (self.leaf_hash.count(), self.compress.count(), self.sponge.count())
    }
}

impl fmt::Display for HashCounts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (l, c, s) = self.snapshot();
        write!(
            f,
            "leaf-hashes: {l}, compressions: {c}, sponge-perms: {s}  (total hash ops: {})",
            l + c + s
        )
    }
}

// ---------------------------------------------------------------------------
// Static reports
// ---------------------------------------------------------------------------

/// Shape of one AIR / table in a batch.
#[derive(Clone, Debug)]
pub struct AirShape {
    pub name: &'static str,
    /// `log2(#rows)`.
    pub log_rows: usize,
    /// Main (witness) trace columns.
    pub main_cols: usize,
    /// Preprocessed (fixed) columns.
    pub preproc_cols: usize,
    /// Number of public values.
    pub num_public: usize,
}

impl AirShape {
    #[must_use]
    pub const fn rows(&self) -> usize {
        1 << self.log_rows
    }

    /// Main-trace cells (`rows × main_cols`).
    #[must_use]
    pub const fn main_cells(&self) -> usize {
        self.rows() * self.main_cols
    }

    /// Preprocessed cells (`rows × preproc_cols`).
    #[must_use]
    pub const fn preproc_cells(&self) -> usize {
        self.rows() * self.preproc_cols
    }
}

/// The circuit shape: one entry per AIR/chip in the batch.
#[derive(Clone, Debug, Default)]
pub struct CircuitShape {
    pub airs: Vec<AirShape>,
}

impl CircuitShape {
    #[must_use]
    pub const fn new(airs: Vec<AirShape>) -> Self {
        Self { airs }
    }

    #[must_use]
    pub fn total_main_cells(&self) -> usize {
        self.airs.iter().map(AirShape::main_cells).sum()
    }

    #[must_use]
    pub fn total_preproc_cells(&self) -> usize {
        self.airs.iter().map(AirShape::preproc_cells).sum()
    }
}

impl fmt::Display for CircuitShape {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "circuit: {} AIR(s)", self.airs.len())?;
        for (i, a) in self.airs.iter().enumerate() {
            writeln!(
                f,
                "  [{i}] {}: {} rows x {} cols (main){}{}",
                a.name,
                a.rows(),
                a.main_cols,
                if a.preproc_cols > 0 {
                    alloc::format!(", {} preprocessed cols", a.preproc_cols)
                } else {
                    alloc::string::String::new()
                },
                if a.num_public > 0 {
                    alloc::format!(", {} public values", a.num_public)
                } else {
                    alloc::string::String::new()
                },
            )?;
        }
        write!(
            f,
            "  totals: {} main cells, {} preprocessed cells",
            self.total_main_cells(),
            self.total_preproc_cells()
        )
    }
}

/// Committed-polynomial + commitment parameters for the tensor/Brakedown PCS,
/// and the two sumcheck round counts.
#[derive(Clone, Debug)]
pub struct PcsInfo {
    /// Variables of the stacked committed multilinear (`log2(total cells)`).
    pub num_variables: usize,
    /// Tensor row split (`rows = 2^log_rows`).
    pub log_rows: usize,
    /// Tensor message width per row (`cols`, a power of two).
    pub cols: usize,
    /// Codeword length `n` per encoded row.
    pub codeword_len: usize,
    /// Column queries opened in the proximity/eval test.
    pub num_queries: usize,
    /// Zerocheck sumcheck rounds (`log2(max #rows)`).
    pub zerocheck_rounds: usize,
    /// Opening-reduction sumcheck rounds (`= num_variables`).
    pub opening_rounds: usize,
}

impl PcsInfo {
    /// Merkle leaves = codeword length (each leaf = one encoded column).
    #[must_use]
    pub const fn num_leaves(&self) -> usize {
        self.codeword_len
    }

    /// Merkle tree depth (`ceil(log2(#leaves))`).
    #[must_use]
    pub fn tree_depth(&self) -> usize {
        let n = self.codeword_len.max(1);
        (usize::BITS - (n - 1).leading_zeros()) as usize
    }
}

impl fmt::Display for PcsInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "poly: {} vars ({} cells); tensor: {} rows x {} cols -> codeword {} ({} leaves, depth {}); {} queries; sumcheck rounds: zerocheck {}, opening {}",
            self.num_variables,
            1usize << self.num_variables,
            1usize << self.log_rows,
            self.cols,
            self.codeword_len,
            self.num_leaves(),
            self.tree_depth(),
            self.num_queries,
            self.zerocheck_rounds,
            self.opening_rounds,
        )
    }
}

/// Serialized byte size of a proof (via `postcard`). `0` if serialization fails.
#[must_use]
pub fn serialized_size<P: Serialize>(proof: &P) -> usize {
    postcard::to_allocvec(proof).map_or(0, |v| v.len())
}

#[cfg(test)]
mod tests {
    use alloc::vec;

    use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
    use p3_symmetric::{PaddingFreeSponge, Permutation, TruncatedPermutation};
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    use super::*;

    fn perm() -> Poseidon2BabyBear<16> {
        Poseidon2BabyBear::new_from_rng_128(&mut SmallRng::seed_from_u64(1))
    }

    #[test]
    fn counting_perm_counts_each_permute() {
        let cp = CountingPerm::new(perm());
        let counter = cp.counter();
        let mut state = [BabyBear::default(); 16];
        for _ in 0..7 {
            cp.permute_mut(&mut state);
        }
        let _ = cp.permute([BabyBear::default(); 16]);
        assert_eq!(counter.count(), 8);
        counter.reset();
        assert_eq!(counter.count(), 0);
    }

    #[test]
    fn counting_compress_counts_each_compression() {
        let compress = TruncatedPermutation::<_, 2, 8, 16>::new(perm());
        let cc = CountingCompress::new(compress);
        let counter = cc.counter();
        let a = [BabyBear::default(); 8];
        for _ in 0..5 {
            let _ = cc.compress([a, a]);
        }
        assert_eq!(counter.count(), 5);
    }

    #[test]
    fn counting_hasher_counts_each_hash() {
        let hasher = PaddingFreeSponge::<_, 16, 8, 8>::new(perm());
        let ch = CountingHasher::new(hasher);
        let counter = ch.counter();
        let input = vec![BabyBear::default(); 20];
        // hash_slice / hash_iter both route through hash_iter -> one count each.
        let _: [BabyBear; 8] = ch.hash_slice(&input);
        let _: [BabyBear; 8] = ch.hash_iter(input.iter().copied());
        assert_eq!(counter.count(), 2);
    }

    #[test]
    fn reports_and_size_render() {
        let shape = CircuitShape::new(vec![AirShape {
            name: "Fib",
            log_rows: 8,
            main_cols: 2,
            preproc_cols: 0,
            num_public: 3,
        }]);
        assert_eq!(shape.total_main_cells(), 256 * 2);
        // Display must not panic and must mention the AIR.
        assert!(alloc::format!("{shape}").contains("Fib"));

        let info = PcsInfo {
            num_variables: 9,
            log_rows: 4,
            cols: 32,
            codeword_len: 55,
            num_queries: 32,
            zerocheck_rounds: 8,
            opening_rounds: 9,
        };
        assert_eq!(info.num_leaves(), 55);
        assert_eq!(info.tree_depth(), 6); // ceil(log2(55)) = 6
        assert!(alloc::format!("{info}").contains("codeword 55"));

        assert!(serialized_size(&vec![1u32, 2, 3]) > 0);
    }
}
