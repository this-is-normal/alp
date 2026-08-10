use core::borrow::BorrowMut;

use p3_field::Field;
use p3_matrix::dense::RowMajorMatrix;

use crate::add::add_points;
use crate::scalarmul::columns::{EcScalarMulCols, NUM_EC_SCALARMUL_COLS};

/// Generates a double-and-add trace for `scalar * base`, one row per bit of `scalar_bits`
/// (least-significant bit first). `scalar_bits.len()` must be a power of two. `base` and
/// `scalar` are pure witness in this AIR (no public values); the final `scalar * base` result
/// ends up in the last row's `result_x`/`result_y` columns.
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

    let mut trace_values = F::zero_vec(num_rows * NUM_EC_SCALARMUL_COLS);

    let mut x_acc = F::ZERO;
    let mut y_acc = F::ONE;
    let (mut x_dbl, mut y_dbl) = base;
    let mut weight = F::ONE;
    let mut scalar_acc = F::ZERO;

    for (row, chunk) in trace_values.chunks_mut(NUM_EC_SCALARMUL_COLS).enumerate() {
        let bit = scalar_bits[row];
        let bit_f = F::from_bool(bit);

        let add_p = x_acc * x_dbl;
        let add_q = y_acc * y_dbl;
        let add_m = d * add_p * add_q;
        let (sum_x, sum_y) = add_points(d, (x_acc, y_acc), (x_dbl, y_dbl));

        let dbl_p = x_dbl * x_dbl;
        let dbl_q = y_dbl * y_dbl;
        let dbl_m = d * dbl_p * dbl_q;
        let (dbl_x3, dbl_y3) = add_points(d, (x_dbl, y_dbl), (x_dbl, y_dbl));

        scalar_acc += bit_f * weight;
        let (result_x, result_y) = if bit { (sum_x, sum_y) } else { (x_acc, y_acc) };

        let cols: &mut EcScalarMulCols<F> = chunk.borrow_mut();
        *cols = EcScalarMulCols {
            bit: bit_f,
            x_acc,
            y_acc,
            x_dbl,
            y_dbl,
            sum_x,
            sum_y,
            add_p,
            add_q,
            add_m,
            dbl_x1: x_dbl,
            dbl_y1: y_dbl,
            dbl_x2: x_dbl,
            dbl_y2: y_dbl,
            dbl_x3,
            dbl_y3,
            dbl_p,
            dbl_q,
            dbl_m,
            weight,
            scalar_acc,
            result_x,
            result_y,
        };

        x_acc = result_x;
        y_acc = result_y;
        x_dbl = dbl_x3;
        y_dbl = dbl_y3;
        weight *= F::TWO;
    }

    RowMajorMatrix::new(trace_values, NUM_EC_SCALARMUL_COLS)
}
