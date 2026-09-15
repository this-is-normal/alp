use core::borrow::BorrowMut;

use alloc::vec;
use alloc::vec::Vec;

use num_bigint::BigUint;
use p3_field::Field;
use p3_matrix::dense::RowMajorMatrix;
use p3_util::log2_strict_usize;

use crate::add::{Point, add_points};
use crate::eddsa::columns::{
    EdDsaPreprocessedCols, EdDsaSharedCols, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS, add_cols,
    add_cols_mut, ladder_cols, ladder_cols_mut, sgc_cols, sgc_cols_mut, shared_cols_mut,
};
use crate::{add, subgroupcheck};

/// Trace height for a `num_bits`-bit scalar: the `num_bits` ladder rows plus the three aux rows
/// (subgroup check on `A`, subgroup check on `R`, the final addition), rounded up to a power of
/// two. Anything above `num_bits + 3` is padding with every selector zero.
///
/// Real ed25519 lands exactly on the boundary: `253 + 3 = 256`.
pub const fn trace_height(num_bits: usize) -> usize {
    (num_bits + 3).next_power_of_two()
}

/// Index of the row carrying the subgroup check on `A`.
pub const fn sgc_a_row(num_bits: usize) -> usize {
    num_bits
}

/// Index of the row carrying the subgroup check on `R`.
pub const fn sgc_r_row(num_bits: usize) -> usize {
    num_bits + 1
}

/// Index of the row carrying `R + [k]A` and the final equality against `[s]B`.
pub const fn add_row(num_bits: usize) -> usize {
    num_bits + 2
}

/// Reconstructs the field element `sum(bit_i * 2^i)` from an LSB-first bit vector, matching how
/// the ladders' `scalar_acc` column accumulates. Deriving `s`/`k` this way (rather than taking
/// them as separate parameters) guarantees they cannot disagree with what the trace reconstructs.
pub fn bits_to_field<F: Field>(bits: &[bool]) -> F {
    let mut acc = F::ZERO;
    let mut weight = F::ONE;
    for &bit in bits {
        if bit {
            acc += weight;
        }
        weight *= F::TWO;
    }
    acc
}

/// LSB-first double-and-add, built on [`add_points`] so it stays valid through the identity
/// (matching the ladders' own convention).
pub fn double_and_add<F: Field>(d: F, base: Point<F>, bits: &[bool]) -> Point<F> {
    let mut acc = (F::ZERO, F::ONE);
    let mut cur = base;
    for &bit in bits {
        if bit {
            acc = add_points(d, acc, cur);
        }
        cur = add_points(d, cur, cur);
    }
    acc
}

/// Panics unless the scalar bit lengths are usable: equal, non-empty, and short enough that
/// `scalar_acc = sum(bit_i * 2^i)` cannot wrap the field characteristic (which would stop the
/// field-equality check against `s`/`k` from uniquely determining the bit vector).
fn check_bit_lengths<F: Field>(s_bits: &[bool], k_bits: &[bool]) -> usize {
    assert_eq!(
        s_bits.len(),
        k_bits.len(),
        "both ladders share the same rows, so `s` and `k` must have the same bit length"
    );
    let num_bits = s_bits.len();
    assert!(num_bits > 0, "scalars must have at least one bit");
    assert!(
        num_bits < F::bits(),
        "scalar bit length ({num_bits}) must be strictly less than the field's bit length ({}); \
         otherwise scalar_acc = sum(bit_i * 2^i) can wrap around the field characteristic, so the \
         field-equality check against `s`/`k` no longer uniquely determines the bit vector",
        F::bits()
    );
    num_bits
}

/// The public values, in the order [`crate::eddsa::air::EdDsaAir`] reads them on the first row:
/// `[a_x, a_y, r_x, r_y, s, k]`. `s` and `k` are derived from the same bits the trace consumes.
pub fn public_values<F: Field>(
    a: Point<F>,
    r: Point<F>,
    s_bits: &[bool],
    k_bits: &[bool],
) -> Vec<F> {
    check_bit_lengths::<F>(s_bits, k_bits);
    vec![
        a.0,
        a.1,
        r.0,
        r.1,
        bits_to_field(s_bits),
        bits_to_field(k_bits),
    ]
}

/// Generates the preprocessed trace for `(d, base_b, num_bits)`: the fixed-base doubling ladder
/// plus the row-type selectors. Depends on no witness, so it is shared by every proof against
/// this base and scalar length.
///
/// `dbl_x, dbl_y` hold `2^i * base_b` on ladder row `i` and the identity elsewhere — junk that no
/// constraint reads, since every constraint touching them is gated by `is_ladder`.
pub fn generate_preprocessed_trace<F: Field>(
    d: F,
    base_b: Point<F>,
    num_bits: usize,
) -> RowMajorMatrix<F> {
    generate_preprocessed_trace_batch(d, base_b, num_bits, 1)
}

/// Stacked ("vertical") version of [`generate_preprocessed_trace`]: tiles the single-instance
/// selector pattern `num_instances` times so `num_instances` signatures share ONE table of
/// `num_instances * trace_height(num_bits)` rows.
///
/// Two things make the tiling correct: the fixed-base doubling chain (`2^i * B`) is **restarted**
/// at the top of every instance, and each instance is bracketed by `is_instance_first`
/// (its row 0, gating the ladder boundary) and `is_instance_last` (its final row, gating OFF the
/// shared-column propagation into the next instance). Every other selector is purely local
/// (`row % h`), so all the row-type constraints tile for free.
pub fn generate_preprocessed_trace_batch<F: Field>(
    d: F,
    base_b: Point<F>,
    num_bits: usize,
    num_instances: usize,
) -> RowMajorMatrix<F> {
    assert!(num_instances >= 1, "need at least one instance");
    let h = trace_height(num_bits);
    let width = NUM_EDDSA_PREPROCESSED_COLS;
    let mut trace_values = F::zero_vec(h * num_instances * width);

    // Build ONE instance, then copy it. Every row depends only on `row % h` —
    // the doubling chain restarts at each instance — so the tail is a memcpy
    // rather than `num_instances` more runs of the fixed-base ladder. That
    // matters because the prover regenerates this per proof when the verifier
    // evaluates the preprocessed trace instead of opening a commitment to it.
    let (mut x_dbl, mut y_dbl) = base_b;
    for (local, chunk) in trace_values[..h * width].chunks_mut(width).enumerate() {
        let cols: &mut EdDsaPreprocessedCols<F> = chunk.borrow_mut();
        let is_ladder = local < num_bits;

        let (dbl_x, dbl_y) = if is_ladder {
            let cur = (x_dbl, y_dbl);
            let (next_x, next_y) = add_points(d, cur, cur);
            x_dbl = next_x;
            y_dbl = next_y;
            cur
        } else {
            (F::ZERO, F::ONE)
        };

        *cols = EdDsaPreprocessedCols {
            dbl_x,
            dbl_y,
            is_ladder: F::from_bool(is_ladder),
            is_ladder_step: F::from_bool(local + 1 < num_bits),
            is_ladder_last: F::from_bool(local + 1 == num_bits),
            is_sgc_a: F::from_bool(local == sgc_a_row(num_bits)),
            is_sgc_r: F::from_bool(local == sgc_r_row(num_bits)),
            is_add: F::from_bool(local == num_bits + 2),
            is_instance_first: F::from_bool(local == 0),
            is_instance_last: F::from_bool(local + 1 == h),
        };
    }
    for i in 1..num_instances {
        let (head, tail) = trace_values.split_at_mut(i * h * width);
        tail[..h * width].copy_from_slice(&head[..h * width]);
    }

    RowMajorMatrix::new(trace_values, width)
}

/// Generates the main trace proving non-cofactored verification of `(R, s)` under public key `A`
/// with challenge scalar `k`, where `s_bits`/`k_bits` are LSB-first and of equal length.
///
/// `l` is the order of the prime-order subgroup, used only to derive the subgroup-check witnesses
/// `(8^-1 mod l) * A` and `(8^-1 mod l) * R`.
///
/// The trace is well-formed for any inputs; it satisfies the AIR only when the signature actually
/// verifies (and `A`, `R` are genuinely in the subgroup). Producing a deliberately invalid trace
/// is how the negative tests work.
pub fn generate_trace_rows<F: Field>(
    d: F,
    base_b: Point<F>,
    l: &BigUint,
    a: Point<F>,
    r: Point<F>,
    s_bits: &[bool],
    k_bits: &[bool],
) -> RowMajorMatrix<F> {
    let num_bits = check_bit_lengths::<F>(s_bits, k_bits);
    let height = trace_height(num_bits);
    let mut trace_values = F::zero_vec(height * NUM_EDDSA_COLS);

    let sb = double_and_add(d, base_b, s_bits);
    let ka = double_and_add(d, a, k_bits);

    // -- ladder rows `0 .. num_bits`: both double-and-add ladders in lockstep --------------------
    //
    // Transcribed from `scalarmul::generate_trace_rows` and `fixed_scalarmul::generate_trace_rows`
    // rather than called: those pad to a power of two, whereas here the ladders occupy exactly
    // `num_bits` rows of a taller table.
    let (mut var_x_acc, mut var_y_acc) = (F::ZERO, F::ONE);
    let (mut var_x_dbl, mut var_y_dbl) = a;
    let (mut fixed_x_acc, mut fixed_y_acc) = (F::ZERO, F::ONE);
    let (mut fixed_x_dbl, mut fixed_y_dbl) = base_b;
    let mut weight = F::ONE;
    let mut var_scalar_acc = F::ZERO;
    let mut fixed_scalar_acc = F::ZERO;

    for (row, chunk) in trace_values
        .chunks_mut(NUM_EDDSA_COLS)
        .take(num_bits)
        .enumerate()
    {
        let cols = ladder_cols_mut(chunk);

        // Variable-base ladder: `[k]A`, with a witnessed doubling block.
        let var_bit = k_bits[row];
        let var_bit_f = F::from_bool(var_bit);
        let var = &mut cols.var;
        var.bit = var_bit_f;
        var.x_acc = var_x_acc;
        var.y_acc = var_y_acc;
        var.x_dbl = var_x_dbl;
        var.y_dbl = var_y_dbl;
        var.add_p = var_x_acc * var_x_dbl;
        var.add_q = var_y_acc * var_y_dbl;
        var.add_m = d * var.add_p * var.add_q;
        (var.sum_x, var.sum_y) = add_points(d, (var_x_acc, var_y_acc), (var_x_dbl, var_y_dbl));
        var.dbl_x1 = var_x_dbl;
        var.dbl_y1 = var_y_dbl;
        var.dbl_x2 = var_x_dbl;
        var.dbl_y2 = var_y_dbl;
        var.dbl_p = var_x_dbl * var_x_dbl;
        var.dbl_q = var_y_dbl * var_y_dbl;
        var.dbl_m = d * var.dbl_p * var.dbl_q;
        (var.dbl_x3, var.dbl_y3) = add_points(d, (var_x_dbl, var_y_dbl), (var_x_dbl, var_y_dbl));
        var_scalar_acc += var_bit_f * weight;
        var.weight = weight;
        var.scalar_acc = var_scalar_acc;
        (var.result_x, var.result_y) = if var_bit {
            (var.sum_x, var.sum_y)
        } else {
            (var_x_acc, var_y_acc)
        };

        // Fixed-base ladder: `[s]B`, reading `2^i * B` from the preprocessed trace instead of
        // witnessing it, hence no doubling block.
        let fixed_bit = s_bits[row];
        let fixed_bit_f = F::from_bool(fixed_bit);
        let fixed = &mut cols.fixed;
        fixed.bit = fixed_bit_f;
        fixed.x_acc = fixed_x_acc;
        fixed.y_acc = fixed_y_acc;
        fixed.add_p = fixed_x_acc * fixed_x_dbl;
        fixed.add_q = fixed_y_acc * fixed_y_dbl;
        fixed.add_m = d * fixed.add_p * fixed.add_q;
        (fixed.sum_x, fixed.sum_y) =
            add_points(d, (fixed_x_acc, fixed_y_acc), (fixed_x_dbl, fixed_y_dbl));
        fixed_scalar_acc += fixed_bit_f * weight;
        fixed.weight = weight;
        fixed.scalar_acc = fixed_scalar_acc;
        (fixed.result_x, fixed.result_y) = if fixed_bit {
            (fixed.sum_x, fixed.sum_y)
        } else {
            (fixed_x_acc, fixed_y_acc)
        };

        (var_x_acc, var_y_acc) = (var.result_x, var.result_y);
        (var_x_dbl, var_y_dbl) = (var.dbl_x3, var.dbl_y3);
        (fixed_x_acc, fixed_y_acc) = (fixed.result_x, fixed.result_y);
        // Advances in lockstep with the preprocessed ladder, which is what actually commits to
        // `2^i * B`.
        (fixed_x_dbl, fixed_y_dbl) =
            add_points(d, (fixed_x_dbl, fixed_y_dbl), (fixed_x_dbl, fixed_y_dbl));
        weight *= F::TWO;
    }

    // -- aux rows: the two subgroup checks, then the final addition ------------------------------
    //
    // Each reuses its standalone generator on a single input and copies the resulting row into
    // this table's region, so the witness layout cannot drift from the AIRs the constraints are
    // transcribed from.
    for (target, row) in [(a, sgc_a_row(num_bits)), (r, sgc_r_row(num_bits))] {
        let witness = subgroupcheck::witness_for_target(d, l, target);
        let generated = subgroupcheck::generate_trace_rows(d, &[witness]);
        let start = row * NUM_EDDSA_COLS;
        *sgc_cols_mut(&mut trace_values[start..start + NUM_EDDSA_COLS]) =
            *sgc_cols(&generated.values);
    }

    {
        let generated = add::generate_trace_rows(d, &[(r, ka)]);
        let start = add_row(num_bits) * NUM_EDDSA_COLS;
        *add_cols_mut(&mut trace_values[start..start + NUM_EDDSA_COLS]) =
            *add_cols(&generated.values);
    }

    // -- shared columns: identical on every row, including the aux and padding rows, since the
    // propagation constraint is a plain `when_transition` that runs to the end of the table ------
    let s = bits_to_field(s_bits);
    let k = bits_to_field(k_bits);
    for chunk in trace_values.chunks_mut(NUM_EDDSA_COLS) {
        *shared_cols_mut(chunk) = EdDsaSharedCols {
            a_x: a.0,
            a_y: a.1,
            r_x: r.0,
            r_y: r.1,
            s,
            k,
            sb_x: sb.0,
            sb_y: sb.1,
            ka_x: ka.0,
            ka_y: ka.1,
        };
    }

    RowMajorMatrix::new(trace_values, NUM_EDDSA_COLS)
}

/// Reads `([s]B, [k]A)` back out of the ladders' last row. Convenience for tests that want to
/// tamper with a result without recomputing it.
pub fn ladder_results<F: Field>(
    trace: &RowMajorMatrix<F>,
    num_bits: usize,
) -> (Point<F>, Point<F>) {
    let start = (num_bits - 1) * NUM_EDDSA_COLS;
    let cols = ladder_cols(&trace.values[start..start + NUM_EDDSA_COLS]);
    (
        (cols.fixed.result_x, cols.fixed.result_y),
        (cols.var.result_x, cols.var.result_y),
    )
}

/// Current- and next-row evaluations of every preprocessed column at `point`,
/// **without** committing to the preprocessed trace.
///
/// The preprocessed trace is a tiling: `generate_preprocessed_trace_batch`
/// repeats one `h = trace_height(num_bits)` row pattern once per instance, and
/// `trace_height` is a power of two, so a row index splits cleanly into
/// `(instance, local)` with no padding. Its multilinear extension therefore
/// collapses — writing `point` as `z_hi ‖ z_lo` with `z_lo` the low
/// `log₂ h` coordinates (MSB-first, so the *trailing* ones):
///
/// ```text
///   MLE_T(z_hi ‖ z_lo) = Σ_i eq(z_hi, i) · Σ_l eq(z_lo, l) · P[l]
///                      = 1 · MLE_P(z_lo)
/// ```
///
/// because `Σ_i eq(z_hi, i) = 1` identically. The same holds for the rotated
/// ("next") view, since a global rotation of a tiling is the tiling of the
/// rotated pattern.
///
/// So the work is `O(h · width)` — **independent of the batch size** — against
/// a committed polynomial that grows with it (2^16 evaluations at 16
/// signatures, 2^22 at 1024). The verifier evaluates public, AIR-fixed data
/// itself instead of trusting an opening, which is strictly the stronger
/// position.
///
/// # Arguments
///
/// `point` is the bound point restricted to this table, `log₂(h · instances)`
/// coordinates, most significant first.
///
/// # Panics
///
/// If `point` is shorter than `log₂ h`.
#[must_use]
pub fn preprocessed_evals_at<F: Field>(
    d: F,
    base_b: Point<F>,
    num_bits: usize,
    point: &[F],
) -> (Vec<F>, Vec<F>) {
    let h = trace_height(num_bits);
    let log_h = log2_strict_usize(h);
    assert!(
        point.len() >= log_h,
        "point has {} coordinates, need at least {log_h}",
        point.len()
    );
    // Only the low coordinates survive the collapse.
    let z_lo = &point[point.len() - log_h..];

    let pattern = generate_preprocessed_trace_batch(d, base_b, num_bits, 1);
    let width = NUM_EDDSA_PREPROCESSED_COLS;

    // eq(z_lo, ·) over the pattern's rows, MSB-first.
    let mut eq = vec![F::ONE];
    for &z in z_lo {
        let mut next = Vec::with_capacity(eq.len() * 2);
        for &e in &eq {
            next.push(e * (F::ONE - z));
            next.push(e * z);
        }
        eq = next;
    }
    debug_assert_eq!(eq.len(), h);

    let mut current = F::zero_vec(width);
    let mut next = F::zero_vec(width);
    for row in 0..h {
        let weight = eq[row];
        let cur_row = &pattern.values[row * width..(row + 1) * width];
        // The rotated pattern: row `h - 1` wraps to row 0, which is also what a
        // global rotation of the tiled trace produces at an instance seam.
        let nxt = (row + 1) % h;
        let nxt_row = &pattern.values[nxt * width..(nxt + 1) * width];
        for c in 0..width {
            current[c] += weight * cur_row[c];
            next[c] += weight * nxt_row[c];
        }
    }
    (current, next)
}
