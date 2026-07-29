use core::borrow::BorrowMut;

use p3_field::Field;
use p3_matrix::dense::RowMajorMatrix;

use crate::add::columns::{EcAddCols, NUM_EC_ADD_COLS};

/// An affine curve point `(x, y)`.
pub type Point<F> = (F, F);

/// Adds two twisted-Edwards points (`a = -1`) using the complete addition law.
///
/// Also valid for doubling, i.e. `p1 == p2`.
pub fn add_points<F: Field>(d: F, p1: Point<F>, p2: Point<F>) -> Point<F> {
    let (x1, y1) = p1;
    let (x2, y2) = p2;
    let p = x1 * x2;
    let q = y1 * y2;
    let m = d * p * q;
    let x3 = (x1 * y2 + y1 * x2) * (F::ONE + m).inverse();
    let y3 = (q + p) * (F::ONE - m).inverse();
    (x3, y3)
}

/// Generates a trace with one row per input pair, filling in the addition result and helper
/// witness columns.
pub fn generate_trace_rows<F: Field>(d: F, pairs: &[(Point<F>, Point<F>)]) -> RowMajorMatrix<F> {
    let num_rows = pairs.len().next_power_of_two().max(1);
    let mut trace = RowMajorMatrix::new(F::zero_vec(num_rows * NUM_EC_ADD_COLS), NUM_EC_ADD_COLS);

    let identity = (F::ZERO, F::ONE);
    for (row, chunk) in trace.values.chunks_mut(NUM_EC_ADD_COLS).enumerate() {
        let cols: &mut EcAddCols<F> = chunk.borrow_mut();
        let (p1, p2) = pairs.get(row).copied().unwrap_or((identity, identity));
        let (x1, y1) = p1;
        let (x2, y2) = p2;
        let p = x1 * x2;
        let q = y1 * y2;
        let m = d * p * q;
        let (x3, y3) = add_points(d, p1, p2);
        *cols = EcAddCols {
            x1,
            y1,
            x2,
            y2,
            x3,
            y3,
            p,
            q,
            m,
        };
    }

    trace
}
