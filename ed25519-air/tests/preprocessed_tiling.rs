//! The preprocessed trace is a *tiling*, so a verifier need not open a
//! commitment to it — it can evaluate it directly, in time independent of the
//! batch size.
//!
//! `generate_preprocessed_trace_batch` tiles one `h = trace_height(num_bits)`
//! row pattern `2^k` times, and `trace_height` is a power of two, so writing a
//! row index as `(instance, local)` gives `T[i·h + l] = P[l]` with no padding.
//! Its multilinear extension therefore collapses:
//!
//! ```text
//!   MLE_T(z_hi ‖ z_lo) = Σ_i eq(z_hi, i) · Σ_l eq(z_lo, l) · P[l]
//!                      = 1 · MLE_P(z_lo)
//! ```
//!
//! because `Σ_i eq(z_hi, i) = 1` identically. The same holds for the rotated
//! ("next") trace, since a global rotation of a tiling is the tiling of the
//! rotated pattern.
//!
//! If these assertions hold, the preprocessed commitment and its opening —
//! roughly half the OpenVM recursion guest, and 42% of the proof blob — are
//! replaceable by `O(h · width)` verifier work.

use p3_ed25519::Ed25519;
use p3_ed25519_air::eddsa::{
    NUM_EDDSA_PREPROCESSED_COLS, generate_preprocessed_trace_batch, trace_height,
};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;

type F = Ed25519;

const NUM_BITS: usize = 253;

fn curve() -> (F, (F, F)) {
    let fe = |s: &str| F::from_biguint(s.parse::<num_bigint::BigUint>().unwrap()).unwrap();
    (
        fe("37095705934669439343138083508754565189542113879843219016388785533085940283555"),
        (
            fe("15112221349535400772501151409588531511454012693041857206046113283949847762202"),
            fe("46316835694926478169428394003475163141307993866256225615783033603165251855960"),
        ),
    )
}

/// `eq(z, ·)` over `z.len()` variables, MSB-first (variable 0 is the most
/// significant bit of the index), matching the fork's convention.
fn eq_table(z: &[F]) -> Vec<F> {
    let mut table = vec![F::ONE];
    for &zi in z {
        let mut next = Vec::with_capacity(table.len() * 2);
        for &t in &table {
            next.push(t * (F::ONE - zi));
            next.push(t * zi);
        }
        table = next;
    }
    table
}

/// Every column's MLE at `z`, by brute force over the whole trace.
fn column_evals(trace: &RowMajorMatrix<F>, z: &[F]) -> Vec<F> {
    let eq = eq_table(z);
    assert_eq!(eq.len(), trace.height());
    (0..trace.width())
        .map(|c| {
            (0..trace.height())
                .map(|r| eq[r] * trace.get(r, c).unwrap())
                .sum()
        })
        .collect()
}

/// The trace rotated by one row, which is what the AIR's `next` view sees.
fn rotate(trace: &RowMajorMatrix<F>) -> RowMajorMatrix<F> {
    let (h, w) = (trace.height(), trace.width());
    let mut values = F::zero_vec(h * w);
    for r in 0..h {
        for c in 0..w {
            values[r * w + c] = trace.get((r + 1) % h, c).unwrap();
        }
    }
    RowMajorMatrix::new(values, w)
}

#[test]
fn trace_height_is_a_power_of_two() {
    // Without this there would be padding rows and the tiling identity would
    // not hold exactly.
    let h = trace_height(NUM_BITS);
    assert!(h.is_power_of_two(), "h = {h}");
    assert_eq!(h, 256);
}

#[test]
fn preprocessed_mle_equals_the_single_instance_pattern() {
    let (d, base_b) = curve();
    let h = trace_height(NUM_BITS);
    let pattern = generate_preprocessed_trace_batch(d, base_b, NUM_BITS, 1);
    assert_eq!(pattern.height(), h);
    assert_eq!(pattern.width(), NUM_EDDSA_PREPROCESSED_COLS);

    let mut rng = ChaCha8Rng::seed_from_u64(0xE0D5A);
    for log_instances in 0..5 {
        let instances = 1 << log_instances;
        let trace = generate_preprocessed_trace_batch(d, base_b, NUM_BITS, instances);
        let num_vars = (h * instances).trailing_zeros() as usize;

        let point: Vec<F> = (0..num_vars).map(|_| rng.random::<F>()).collect();
        // MSB-first: the leading variables select the instance, the trailing 8
        // the row within it.
        let z_lo = &point[num_vars - 8..];

        assert_eq!(
            column_evals(&trace, &point),
            column_evals(&pattern, z_lo),
            "current-row evaluations disagree at {instances} instances"
        );
        assert_eq!(
            column_evals(&rotate(&trace), &point),
            column_evals(&rotate(&pattern), z_lo),
            "next-row evaluations disagree at {instances} instances"
        );
    }
}

/// `preprocessed_evals_at` must return exactly what a commitment opening would:
/// every column's MLE over the **full** trace, current and next, at the bound
/// point. This is the identity the whole optimisation rests on.
#[test]
fn direct_evaluation_matches_the_committed_trace() {
    use p3_ed25519_air::eddsa::preprocessed_evals_at;

    let (d, base_b) = curve();
    let h = trace_height(NUM_BITS);
    let mut rng = ChaCha8Rng::seed_from_u64(0x1DEA);

    for log_instances in 0..5 {
        let instances = 1 << log_instances;
        let trace = generate_preprocessed_trace_batch(d, base_b, NUM_BITS, instances);
        let num_vars = (h * instances).trailing_zeros() as usize;
        let point: Vec<F> = (0..num_vars).map(|_| rng.random::<F>()).collect();

        let (current, next) = preprocessed_evals_at(d, base_b, NUM_BITS, &point);

        assert_eq!(
            current,
            column_evals(&trace, &point),
            "current-row evaluations wrong at {instances} instances"
        );
        assert_eq!(
            next,
            column_evals(&rotate(&trace), &point),
            "next-row evaluations wrong at {instances} instances"
        );
    }
}
