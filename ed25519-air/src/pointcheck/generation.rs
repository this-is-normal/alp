use core::borrow::BorrowMut;

use p3_field::Field;
use p3_matrix::dense::RowMajorMatrix;

use crate::pointcheck::columns::{EcPointCheckCols, NUM_EC_POINTCHECK_COLS};

/// Generates a trace with one row per input point, filling in the helper witness columns.
///
/// Pads with the identity point `(0, 1)`, which is on-curve for any `a`, `d`.
pub fn generate_trace_rows<F: Field>(points: &[(F, F)]) -> RowMajorMatrix<F> {
    let num_rows = points.len().next_power_of_two().max(1);
    let mut trace = RowMajorMatrix::new(
        F::zero_vec(num_rows * NUM_EC_POINTCHECK_COLS),
        NUM_EC_POINTCHECK_COLS,
    );

    let identity = (F::ZERO, F::ONE);
    for (row, chunk) in trace.values.chunks_mut(NUM_EC_POINTCHECK_COLS).enumerate() {
        let cols: &mut EcPointCheckCols<F> = chunk.borrow_mut();
        let (x, y) = points.get(row).copied().unwrap_or(identity);
        let x2 = x * x;
        let y2 = y * y;
        *cols = EcPointCheckCols { x, y, x2, y2 };
    }

    trace
}
