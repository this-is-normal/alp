use core::borrow::Borrow;

use p3_air::{Air, AirBuilder, BaseAir, WindowAccess};
use p3_field::{Dup, PrimeCharacteristicRing};

use crate::pointcheck::columns::{EcPointCheckCols, NUM_EC_POINTCHECK_COLS};

/// AIR checking that each row's `(x, y)` lies on the twisted-Edwards curve
/// `-x^2 + y^2 = 1 + d*x^2*y^2` (`a = -1`, matching ed25519 and [`EcAddAir`][crate::add::EcAddAir]).
///
/// Each row is an independent check; there are no transition or boundary constraints.
#[derive(Debug, Clone, Copy)]
pub struct EcPointCheckAir<F> {
    pub d: F,
}

impl<F: Sync> BaseAir<F> for EcPointCheckAir<F> {
    fn width(&self) -> usize {
        NUM_EC_POINTCHECK_COLS
    }
}

impl<AB: AirBuilder> Air<AB> for EcPointCheckAir<AB::F> {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local: &EcPointCheckCols<AB::Var> = main.current_slice().borrow();
        let &EcPointCheckCols { x, y, x2, y2 } = local;

        // x2 = x*x, y2 = y*y
        builder.assert_zero(x * x - x2);
        builder.assert_zero(y * y - y2);

        // -x2 + y2 == 1 + d*x2*y2
        let lhs = y2 - x2;
        let rhs = AB::Expr::ONE + (x2 * y2) * self.d.dup();
        builder.assert_zero(lhs - rhs);
    }
}
