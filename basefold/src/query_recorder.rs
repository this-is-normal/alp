//! Record which foldable-code indices a verifier reads.
//!
//! **Instrumentation, not a verifier feature.** A recursion guest needs
//! `twist[idx]` and `1/twist[idx]` at the `O(q·log N)` positions its queries
//! land on, but it cannot know which those are until it runs — and the host
//! can, because verification is deterministic. Recording them here lets the
//! host ship exactly those elements, rather than the guest carrying the whole
//! table in its program image (16.8 MB of table is ~86 MB of `.vmexe`, and
//! image size is charged on every `prove app`).
//!
//! Enabled by the `record-queries` feature, which pulls in `std` for the
//! thread-local. It is never enabled in a guest build.

#[cfg(feature = "record-queries")]
mod imp {
    extern crate std;

    use alloc::vec::Vec;
    use std::cell::RefCell;

    std::thread_local! {
        static QUERIES: RefCell<Vec<usize>> = const { RefCell::new(Vec::new()) };
    }

    /// Note that `idx` was read. Called from the code's lazy accessors.
    pub fn record(idx: usize) {
        QUERIES.with(|q| q.borrow_mut().push(idx));
    }

    /// Drop everything recorded so far.
    pub fn reset() {
        QUERIES.with(|q| q.borrow_mut().clear());
    }

    /// The recorded indices, **sorted and deduplicated** — the order a verifier
    /// happens to touch them in is an implementation detail, and a host that
    /// shipped them in that order would couple the two sides to it. Sorted
    /// means the guest can binary-search by index instead.
    #[must_use]
    pub fn take_sorted() -> Vec<usize> {
        QUERIES.with(|q| {
            let mut v = core::mem::take(&mut *q.borrow_mut());
            v.sort_unstable();
            v.dedup();
            v
        })
    }

    /// Total reads including repeats — the count that predicts guest cost,
    /// where `take_sorted().len()` predicts blob size.
    #[must_use]
    pub fn read_count() -> usize {
        QUERIES.with(|q| q.borrow().len())
    }
}

#[cfg(not(feature = "record-queries"))]
mod imp {
    use alloc::vec::Vec;

    #[inline(always)]
    pub fn record(_idx: usize) {}
    pub fn reset() {}
    #[must_use]
    pub fn take_sorted() -> Vec<usize> {
        Vec::new()
    }
    #[must_use]
    pub fn read_count() -> usize {
        0
    }
}

pub use imp::{read_count, record, reset, take_sorted};
