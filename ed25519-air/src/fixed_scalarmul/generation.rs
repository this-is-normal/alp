use core::borrow::BorrowMut;

use p3_field::Field;
use p3_matrix::dense::RowMajorMatrix;

use crate::add::add_points;
use crate::fixed_scalarmul::columns::{
    EcFixedScalarMulCols, EcFixedScalarMulPreprocessedCols, NUM_EC_FIXED_SCALARMUL_COLS,
    NUM_EC_FIXED_SCALARMUL_PREPROCESSED_COLS,
};

/// Generates the preprocessed doubling ladder `2^i * base` for `i` in `0..num_rows`, one
/// `(x, y)` row each. Depends only on `(base, num_rows)`, not on any witness — shared across
/// every proof against this base.
pub fn generate_preprocessed_trace<F: Field>(
    d: F,
    base: (F, F),
    num_rows: usize,
) -> RowMajorMatrix<F> {
    assert!(
        num_rows.is_power_of_two(),
        "preprocessed trace height must be a power of two"
    );

    let mut trace_values = F::zero_vec(num_rows * NUM_EC_FIXED_SCALARMUL_PREPROCESSED_COLS);

    let (mut x_dbl, mut y_dbl) = base;
    for chunk in trace_values.chunks_mut(NUM_EC_FIXED_SCALARMUL_PREPROCESSED_COLS) {
        let cols: &mut EcFixedScalarMulPreprocessedCols<F> = chunk.borrow_mut();
        *cols = EcFixedScalarMulPreprocessedCols { x_dbl, y_dbl };

        let (next_x, next_y) = add_points(d, (x_dbl, y_dbl), (x_dbl, y_dbl));
        x_dbl = next_x;
        y_dbl = next_y;
    }

    RowMajorMatrix::new(trace_values, NUM_EC_FIXED_SCALARMUL_PREPROCESSED_COLS)
}

/// Generates a double-and-add trace for `scalar * base`, one row per bit of `scalar_bits`
/// (least-significant bit first). `scalar_bits.len()` must be a power of two and must match
/// the `num_rows` used to build the corresponding preprocessed trace. `scalar` is pure witness
/// in this AIR (no public values); the final `scalar * base` result ends up in the last row's
/// `result_x`/`result_y` columns.
pub fn generate_trace_rows<F: Field>(
    d: F,
    base: (F, F),
    scalar_bits: &[bool],
) -> RowMajorMatrix<F> {
    let num_rows = scalar_bits.len();
    assert!(
        num_rows.is_power_of_two(),
        "trace height must be a power of two"
    );
    assert!(
        num_rows < F::bits(),
        "scalar_bits.len() ({num_rows}) must be strictly less than the field's bit length \
         ({}); otherwise scalar_acc = sum(bit_i * 2^i) can wrap around the field \
         characteristic, so the field-equality check against `scalar` no longer uniquely \
         determines the bit vector",
        F::bits()
    );

    let mut trace_values = F::zero_vec(num_rows * NUM_EC_FIXED_SCALARMUL_COLS);

    let mut x_acc = F::ZERO;
    let mut y_acc = F::ONE;
    let (mut x_dbl, mut y_dbl) = base;
    let mut weight = F::ONE;
    let mut scalar_acc = F::ZERO;

    for (row, chunk) in trace_values
        .chunks_mut(NUM_EC_FIXED_SCALARMUL_COLS)
        .enumerate()
    {
        let bit = scalar_bits[row];
        let bit_f = F::from_bool(bit);

        let add_p = x_acc * x_dbl;
        let add_q = y_acc * y_dbl;
        let add_m = d * add_p * add_q;
        let (sum_x, sum_y) = add_points(d, (x_acc, y_acc), (x_dbl, y_dbl));

        scalar_acc += bit_f * weight;
        let (result_x, result_y) = if bit { (sum_x, sum_y) } else { (x_acc, y_acc) };

        let cols: &mut EcFixedScalarMulCols<F> = chunk.borrow_mut();
        *cols = EcFixedScalarMulCols {
            bit: bit_f,
            x_acc,
            y_acc,
            sum_x,
            sum_y,
            add_p,
            add_q,
            add_m,
            weight,
            scalar_acc,
            result_x,
            result_y,
        };

        x_acc = result_x;
        y_acc = result_y;
        // Advances in lockstep with `generate_preprocessed_trace`'s ladder; not itself
        // witnessed, since the preprocessed trace already commits to `2^i * base`.
        let (next_x, next_y) = add_points(d, (x_dbl, y_dbl), (x_dbl, y_dbl));
        x_dbl = next_x;
        y_dbl = next_y;
        weight *= F::TWO;
    }

    RowMajorMatrix::new(trace_values, NUM_EC_FIXED_SCALARMUL_COLS)
}
