//! The recorder must see exactly the indices a lazy verifier would need.
//!
//! This is what lets a host ship those elements instead of a guest baking the
//! whole table. If it under-reports, the guest starves; if the indices are
//! wrong, verification fails with a transcript mismatch that looks like a bug
//! anywhere else.

#![cfg(feature = "record-queries")]

use p3_basefold::{FoldableCode, RandomFoldableCode, level_base, query_recorder, sample_nonzero_at};
use p3_ed25519::Ed25519;

type F = Ed25519;

const SEED: u64 = 0xBA5E_F01D;
const RATE: usize = 4;
const LOG_LEN: usize = 6;

/// Reading a position through a **materialised** code records the same index a
/// lazy one would derive. This is the property that lets the host record using
/// its ordinary prover-side config, with no separate lazy verifier build.
#[test]
fn a_materialised_code_records_what_a_lazy_one_would_derive() {
    let eager = RandomFoldableCode::<F>::new(LOG_LEN, RATE, SEED);
    let lazy = RandomFoldableCode::<F>::new_verifier(LOG_LEN, RATE, SEED);

    for (level, j) in [(1usize, 0usize), (3, 5), (LOG_LEN, 17)] {
        query_recorder::reset();
        let a = eager.diag_at(level, j);
        let from_eager = query_recorder::take_sorted();

        query_recorder::reset();
        let b = lazy.diag_at(level, j);
        let from_lazy = query_recorder::take_sorted();

        assert_eq!(a, b, "eager and lazy must agree on the value");
        assert_eq!(from_eager, from_lazy, "and on the recorded index");
        assert_eq!(from_eager, vec![level_base(RATE, level) + j]);
    }
}

/// Every recorded index must resolve, through `sample_nonzero_at`, to the value
/// the code actually uses — that equivalence is what makes a shipped element or
/// a baked table substitutable for derivation.
#[test]
fn recorded_indices_resolve_to_the_codes_own_values() {
    let code = RandomFoldableCode::<F>::new(LOG_LEN, RATE, SEED);
    query_recorder::reset();

    for level in 1..=LOG_LEN {
        let n = RATE << (level - 1);
        for j in (0..n).step_by(7) {
            let t = code.diag_at(level, j);
            assert_eq!(
                t,
                sample_nonzero_at::<F>(SEED, level_base(RATE, level) + j),
                "index layout must match sample_nonzero_at at (level {level}, j {j})"
            );
        }
    }
    let idxs = query_recorder::take_sorted();
    assert!(!idxs.is_empty(), "recorder saw nothing");
    assert!(idxs.windows(2).all(|w| w[0] < w[1]), "must be sorted and unique");
}

/// Repeated reads of one position collapse to a single shipped element, but the
/// *read* count is what predicts guest cost. Both numbers matter, and they differ.
#[test]
fn reads_and_distinct_indices_are_counted_separately() {
    let code = RandomFoldableCode::<F>::new(LOG_LEN, RATE, SEED);
    query_recorder::reset();
    for _ in 0..5 {
        let _ = code.diag_at(2, 3);
    }
    assert_eq!(query_recorder::read_count(), 5, "five reads");
    assert_eq!(query_recorder::take_sorted().len(), 1, "one distinct index");
}
