use core::borrow::Borrow;

use p3_air::{Air, AirBuilder, BaseAir, WindowAccess};
use p3_field::{Dup, PrimeCharacteristicRing};

use crate::add::columns::{EcAddCols, NUM_EC_ADD_COLS};

/// AIR for the complete twisted-Edwards addition law (curve parameter `a = -1` fixed,
/// `d` supplied at construction so callers can plug in the real ed25519 base-field value
/// once it lands, or a placeholder for structural testing over another native field).
///
/// Each row is an independent addition; there are no transition or boundary constraints.
#[derive(Debug, Clone, Copy)]
pub struct EcAddAir<F> {
    pub d: F,
}

impl<F: Sync> BaseAir<F> for EcAddAir<F> {
    fn width(&self) -> usize {
        NUM_EC_ADD_COLS
    }
}

impl<AB: AirBuilder> Air<AB> for EcAddAir<AB::F> {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local: &EcAddCols<AB::Var> = main.current_slice().borrow();
        let &EcAddCols {
            x1,
            y1,
            x2,
            y2,
            x3,
            y3,
            p,
            q,
            m,
        } = local;

        // p = x1*x2, q = y1*y2
        builder.assert_zero(x1 * x2 - p);
        builder.assert_zero(y1 * y2 - q);

        // m = d*p*q
        builder.assert_zero((p * q) * self.d.dup() - m);

        // x3*(1+m) = x1*y2 + y1*x2
        let one_plus_m = AB::Expr::ONE + m;
        builder.assert_zero(x3 * one_plus_m - (x1 * y2 + y1 * x2));

        // y3*(1-m) = q + p   (using a = -1, so y1*y2 - a*x1*x2 = q + p)
        let one_minus_m = AB::Expr::ONE - m;
        builder.assert_zero(y3 * one_minus_m - (q + p));
    }
}
