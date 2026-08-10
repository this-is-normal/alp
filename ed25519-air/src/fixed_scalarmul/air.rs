use core::borrow::Borrow;

use p3_air::{Air, AirBuilder, BaseAir, WindowAccess};
use p3_field::{Dup, Field, PrimeCharacteristicRing};
use p3_matrix::dense::RowMajorMatrix;

use crate::fixed_scalarmul::columns::{
    EcFixedScalarMulCols, EcFixedScalarMulPreprocessedCols, NUM_EC_FIXED_SCALARMUL_COLS,
    NUM_EC_FIXED_SCALARMUL_PREPROCESSED_COLS,
};
use crate::fixed_scalarmul::generation::generate_preprocessed_trace;

/// AIR for `scalar * base` via double-and-add over one scalar bit per row (LSB first),
/// against a `base` fixed at circuit-definition time.
///
/// Unlike [`crate::scalarmul::EcScalarMulAir`], `2^i * base` is not re-derived and proved
/// on-circuit every row; it's supplied as a preprocessed column (see
/// [`Self::preprocessed_trace`]), so the AIR embeds the twisted-Edwards add law directly
/// (rather than as an [`crate::add::EcAddAir`] sub-AIR — see `docs-ai/SubAirBuilder.md` for
/// why) exactly once per row, for `acc + dbl` only.
#[derive(Debug, Clone, Copy)]
pub struct EcFixedScalarMulAir<F> {
    pub d: F,
    pub base: (F, F),
    /// Number of rows in the trace this AIR will be evaluated against; needed up front since
    /// `preprocessed_trace` derives the doubling ladder from `(base, num_rows)`.
    pub num_rows: usize,
}

impl<F: Field + Sync> BaseAir<F> for EcFixedScalarMulAir<F> {
    fn width(&self) -> usize {
        NUM_EC_FIXED_SCALARMUL_COLS
    }

    fn preprocessed_trace(&self) -> Option<RowMajorMatrix<F>> {
        Some(generate_preprocessed_trace(
            self.d,
            self.base,
            self.num_rows,
        ))
    }

    fn preprocessed_width(&self) -> usize {
        NUM_EC_FIXED_SCALARMUL_PREPROCESSED_COLS
    }
}

impl<AB: AirBuilder> Air<AB> for EcFixedScalarMulAir<AB::F>
where
    AB::F: Field,
{
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local: &EcFixedScalarMulCols<AB::Var> = main.current_slice().borrow();
        let next: &EcFixedScalarMulCols<AB::Var> = main.next_slice().borrow();
        let prep: &EcFixedScalarMulPreprocessedCols<AB::Var> =
            builder.preprocessed().current_slice().borrow();

        let &EcFixedScalarMulCols {
            bit,
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
        } = local;
        let &EcFixedScalarMulPreprocessedCols { x_dbl, y_dbl } = prep;

        // Add law for `acc + dbl` (mirrors `EcAddAir::eval`), with `x2, y2` read from the
        // preprocessed row instead of a witness column.
        builder.assert_zero(x_acc * x_dbl - add_p);
        builder.assert_zero(y_acc * y_dbl - add_q);
        builder.assert_zero((add_p * add_q) * self.d.dup() - add_m);
        let one_plus_m = AB::Expr::ONE + add_m;
        builder.assert_zero(sum_x * one_plus_m - (x_acc * y_dbl + y_acc * x_dbl));
        let one_minus_m = AB::Expr::ONE - add_m;
        builder.assert_zero(sum_y * one_minus_m - (add_q + add_p));

        // The scalar bit must be boolean.
        builder.assert_bool(bit);

        // First row: accumulator starts at the identity, and the weight starts at 2^0 = 1.
        // (`x_dbl_0 = base` is true by construction of the preprocessed trace, not asserted
        // here.)
        let mut first_row = builder.when_first_row();
        first_row.assert_zero(x_acc);
        first_row.assert_one(y_acc);
        first_row.assert_one(weight);
        first_row.assert_eq(scalar_acc, bit * weight);

        // Transition: fold this row's bit into the accumulator (if set), and advance the
        // weight/running scalar sum. (Advancing the doubling point is likewise true by
        // construction of the preprocessed trace, not constrained here.)
        let mut transition = builder.when_transition();
        transition.assert_eq(next.x_acc, bit * sum_x + (AB::Expr::ONE - bit) * x_acc);
        transition.assert_eq(next.y_acc, bit * sum_y + (AB::Expr::ONE - bit) * y_acc);
        transition.assert_eq(next.weight, weight * AB::Expr::TWO);
        transition.assert_eq(next.scalar_acc, scalar_acc + next.bit * next.weight);

        // Every row: `result_x`/`result_y` hold the accumulator after folding in this row's
        // bit — redundant with the transition constraint above except on the last row, where
        // it's the only witnessed copy of the final `scalar * base` result. `scalar_acc`'s
        // final value plays the same role for the reconstructed scalar. Neither is tied to a
        // public value here; a composing AIR does that itself if it needs to.
        builder.assert_eq(result_x, bit * sum_x + (AB::Expr::ONE - bit) * x_acc);
        builder.assert_eq(result_y, bit * sum_y + (AB::Expr::ONE - bit) * y_acc);
    }
}
