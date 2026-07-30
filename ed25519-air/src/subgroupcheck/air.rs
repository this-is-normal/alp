use core::borrow::Borrow;

use p3_air::{Air, AirBuilder, BaseAir, WindowAccess};
use p3_field::{Dup, Field, PrimeCharacteristicRing};
use p3_uni_stark::SubAirBuilder;

use crate::add::EcAddAir;
use crate::pointcheck::EcPointCheckAir;
use crate::subgroupcheck::columns::{
    CC_BLOCK_START, D1_BLOCK_START, D2_BLOCK_START, D3_BLOCK_START, EcSubgroupCheckCols,
    NUM_EC_SUBGROUPCHECK_COLS,
};

/// AIR proving that `(q_x, q_y)` is `8*P` for some on-curve witness point `P` — a
/// subgroup-membership certificate for the prime-order subgroup of ed25519 (see
/// `docs/plans/SubgroupCheck.md`).
///
/// Fixes `a = -1`, matching [`EcAddAir`], which this AIR embeds three times (one per
/// doubling) alongside one embedding of [`EcPointCheckAir`] (to check `P` itself is on-curve).
#[derive(Debug, Clone, Copy)]
pub struct EcSubgroupCheckAir<F> {
    pub d: F,
}

impl<F: Sync> BaseAir<F> for EcSubgroupCheckAir<F> {
    fn width(&self) -> usize {
        NUM_EC_SUBGROUPCHECK_COLS
    }
}

impl<AB: AirBuilder> Air<AB> for EcSubgroupCheckAir<AB::F>
where
    AB::F: Field,
{
    fn eval(&self, builder: &mut AB) {
        let point_check_air = EcPointCheckAir {
            a: -AB::F::ONE,
            d: self.d.dup(),
        };
        let add_air = EcAddAir { d: self.d.dup() };

        {
            let mut sub_builder = SubAirBuilder::<AB, EcPointCheckAir<AB::F>, AB::Var>::new(
                builder,
                CC_BLOCK_START..CC_BLOCK_START + 4,
            );
            point_check_air.eval(&mut sub_builder);
        }
        {
            let mut sub_builder = SubAirBuilder::<AB, EcAddAir<AB::F>, AB::Var>::new(
                builder,
                D1_BLOCK_START..D1_BLOCK_START + 9,
            );
            add_air.eval(&mut sub_builder);
        }
        {
            let mut sub_builder = SubAirBuilder::<AB, EcAddAir<AB::F>, AB::Var>::new(
                builder,
                D2_BLOCK_START..D2_BLOCK_START + 9,
            );
            add_air.eval(&mut sub_builder);
        }
        {
            let mut sub_builder = SubAirBuilder::<AB, EcAddAir<AB::F>, AB::Var>::new(
                builder,
                D3_BLOCK_START..D3_BLOCK_START + 9,
            );
            add_air.eval(&mut sub_builder);
        }

        let main = builder.main();
        let local: &EcSubgroupCheckCols<AB::Var> = main.current_slice().borrow();
        let &EcSubgroupCheckCols {
            p_x,
            p_y,
            cc_x,
            cc_y,
            d1_x1,
            d1_y1,
            d1_x2,
            d1_y2,
            d1_x3,
            d1_y3,
            d2_x1,
            d2_y1,
            d2_x2,
            d2_y2,
            d2_x3,
            d2_y3,
            d3_x1,
            d3_y1,
            d3_x2,
            d3_y2,
            d3_x3,
            d3_y3,
            q_x,
            q_y,
            ..
        } = local;

        // Feed P into the point-check sub-air block.
        builder.assert_zero(cc_x - p_x);
        builder.assert_zero(cc_y - p_y);

        // First doubling block: both logical inputs are P itself.
        builder.assert_zero(d1_x1 - p_x);
        builder.assert_zero(d1_y1 - p_y);
        builder.assert_zero(d1_x2 - p_x);
        builder.assert_zero(d1_y2 - p_y);

        // Second doubling block: both logical inputs are 2P (the first block's output).
        builder.assert_zero(d2_x1 - d1_x3);
        builder.assert_zero(d2_y1 - d1_y3);
        builder.assert_zero(d2_x2 - d1_x3);
        builder.assert_zero(d2_y2 - d1_y3);

        // Third doubling block: both logical inputs are 4P (the second block's output).
        builder.assert_zero(d3_x1 - d2_x3);
        builder.assert_zero(d3_y1 - d2_y3);
        builder.assert_zero(d3_x2 - d2_x3);
        builder.assert_zero(d3_y2 - d2_y3);

        // The point under test must equal 8P (the third block's output).
        builder.assert_zero(q_x - d3_x3);
        builder.assert_zero(q_y - d3_y3);
    }
}
