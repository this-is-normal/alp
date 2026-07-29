use core::borrow::Borrow;

use p3_air::{Air, AirBuilder, BaseAir, WindowAccess};
use p3_field::{Dup, PrimeCharacteristicRing};
use p3_uni_stark::SubAirBuilder;

use crate::add::EcAddAir;
use crate::scalarmul::columns::{
    ADD_BLOCK_START, DBL_BLOCK_START, EcScalarMulCols, NUM_EC_SCALARMUL_COLS,
};

/// Public values: `[base_x, base_y, scalar, result_x, result_y]`.
pub const NUM_PUBLIC_VALUES: usize = 5;

/// AIR for `scalar * base` via double-and-add over one scalar bit per row (LSB first),
/// built by embedding [`EcAddAir`] as a sub-AIR twice per row (once for the conditional
/// accumulator update, once for doubling).
#[derive(Debug, Clone, Copy)]
pub struct EcScalarMulAir<F> {
    pub d: F,
}

impl<F: Sync> BaseAir<F> for EcScalarMulAir<F> {
    fn width(&self) -> usize {
        NUM_EC_SCALARMUL_COLS
    }

    fn num_public_values(&self) -> usize {
        NUM_PUBLIC_VALUES
    }
}

impl<AB: AirBuilder> Air<AB> for EcScalarMulAir<AB::F> {
    fn eval(&self, builder: &mut AB) {
        let add_air = EcAddAir { d: self.d.dup() };

        {
            let mut sub_builder = SubAirBuilder::<AB, EcAddAir<AB::F>, AB::Var>::new(
                builder,
                ADD_BLOCK_START..ADD_BLOCK_START + 9,
            );
            add_air.eval(&mut sub_builder);
        }
        {
            let mut sub_builder = SubAirBuilder::<AB, EcAddAir<AB::F>, AB::Var>::new(
                builder,
                DBL_BLOCK_START..DBL_BLOCK_START + 9,
            );
            add_air.eval(&mut sub_builder);
        }

        let main = builder.main();
        let local: &EcScalarMulCols<AB::Var> = main.current_slice().borrow();
        let next: &EcScalarMulCols<AB::Var> = main.next_slice().borrow();
        let &EcScalarMulCols {
            bit,
            x_acc,
            y_acc,
            x_dbl,
            y_dbl,
            sum_x,
            sum_y,
            dbl_x1,
            dbl_y1,
            dbl_x2,
            dbl_y2,
            dbl_x3,
            dbl_y3,
            weight,
            scalar_acc,
            ..
        } = local;

        // The doubling sub-air block doubles `(x_dbl, y_dbl)` with itself; wire its two
        // logical inputs to the actual doubling-point columns.
        builder.assert_zero(dbl_x1 - x_dbl);
        builder.assert_zero(dbl_y1 - y_dbl);
        builder.assert_zero(dbl_x2 - x_dbl);
        builder.assert_zero(dbl_y2 - y_dbl);

        // The scalar bit must be boolean.
        builder.assert_bool(bit);

        let public_values = builder.public_values();
        let base_x = public_values[0];
        let base_y = public_values[1];
        let scalar = public_values[2];
        let result_x = public_values[3];
        let result_y = public_values[4];

        // First row: accumulator starts at the identity, doubling point starts at the base
        // point, and the weight starts at 2^0 = 1.
        let mut first_row = builder.when_first_row();
        first_row.assert_zero(x_acc);
        first_row.assert_one(y_acc);
        first_row.assert_eq(x_dbl, base_x);
        first_row.assert_eq(y_dbl, base_y);
        first_row.assert_one(weight);
        first_row.assert_eq(scalar_acc, bit * weight);

        // Transition: fold this row's bit into the accumulator (if set), double the doubling
        // point, and advance the weight/running scalar sum.
        let mut transition = builder.when_transition();
        transition.assert_eq(next.x_acc, bit * sum_x + (AB::Expr::ONE - bit) * x_acc);
        transition.assert_eq(next.y_acc, bit * sum_y + (AB::Expr::ONE - bit) * y_acc);
        transition.assert_eq(next.x_dbl, dbl_x3);
        transition.assert_eq(next.y_dbl, dbl_y3);
        transition.assert_eq(next.weight, weight * AB::Expr::TWO);
        transition.assert_eq(
            next.scalar_acc,
            scalar_acc + next.bit * next.weight,
        );

        // Last row: the scalar bits must reconstruct the public scalar, and folding in the
        // last bit must produce the public result point.
        let mut last_row = builder.when_last_row();
        last_row.assert_eq(scalar_acc, scalar);
        last_row.assert_eq(
            result_x,
            bit * sum_x + (AB::Expr::ONE - bit) * x_acc,
        );
        last_row.assert_eq(
            result_y,
            bit * sum_y + (AB::Expr::ONE - bit) * y_acc,
        );
    }
}
