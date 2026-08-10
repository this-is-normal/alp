use core::borrow::Borrow;

use p3_air::{Air, AirBuilder, BaseAir, WindowAccess};
use p3_field::{Dup, Field, PrimeCharacteristicRing};
use p3_matrix::dense::RowMajorMatrix;

use crate::add::Point;
use crate::eddsa::columns::{
    EdDsaPreprocessedCols, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS, NUM_EDDSA_PUBLIC_VALUES,
    add_cols, ladder_cols, sgc_cols, shared_cols,
};
use crate::eddsa::generation::generate_preprocessed_trace;

/// Monolithic AIR for non-cofactored RFC 8032 verification of `(R, s)` under public key `A` with
/// challenge scalar `k`, over public values `[a_x, a_y, r_x, r_y, s, k]`.
///
/// One table proves all of it:
/// - `[s]B` (fixed-base ladder) and `[k]A` (variable-base ladder) on rows `0 .. num_bits`,
/// - `A` is in the prime-order subgroup on row `num_bits`,
/// - `R` is in the prime-order subgroup on row `num_bits + 1`,
/// - `R + [k]A == [s]B` on row `num_bits + 2`,
///
/// with every "same value" tie between those pieces enforced in-circuit via the shared columns
/// (see [`crate::eddsa::columns::EdDsaSharedCols`]). Rows past `num_bits + 2` are padding.
///
/// The row type is selected by a preprocessed selector, and every constraint below is emitted
/// through `builder.when(selector)`, so a row type's constraints vanish everywhere else — on
/// padding rows all selectors are zero and nothing is enforced but the shared-column propagation.
/// The constraint logic is transcribed from the five standalone AIRs (`EcAddAir`,
/// `EcPointCheckAir`, `EcSubgroupCheckAir`, `EcScalarMulAir`, `EcFixedScalarMulAir`) rather than
/// composed from them: `SubAirBuilder` restricts a sub-AIR to a column range but applies no
/// selector, so the gated versions have to be written out by hand.
///
/// Max constraint degree is 3: the building blocks' constraints are degree 2, and one selector
/// factor makes them degree 3.
///
/// **What is not proved:** that `k == H(dom2 || R || A || M) mod L`, and that `s`/`k` are reduced
/// mod `L`. See "What is proved" in `docs/EdDSA.md`.
#[derive(Debug, Clone, Copy)]
pub struct EdDsaAir<F> {
    /// Twisted-Edwards curve parameter (`a = -1` is fixed).
    pub d: F,
    /// The generator `B`, fixed at circuit-definition time; the fixed-base ladder's doubling
    /// chain `2^i * B` is baked into the preprocessed trace.
    pub base_b: Point<F>,
    /// Scalar bit length `N`; the trace is `next_power_of_two(N + 3)` rows tall.
    pub num_bits: usize,
}

impl<F: Field + Sync> BaseAir<F> for EdDsaAir<F> {
    fn width(&self) -> usize {
        NUM_EDDSA_COLS
    }

    fn preprocessed_trace(&self) -> Option<RowMajorMatrix<F>> {
        Some(generate_preprocessed_trace(
            self.d,
            self.base_b,
            self.num_bits,
        ))
    }

    fn preprocessed_width(&self) -> usize {
        NUM_EDDSA_PREPROCESSED_COLS
    }

    fn num_public_values(&self) -> usize {
        NUM_EDDSA_PUBLIC_VALUES
    }
}

/// Emits the complete twisted-Edwards addition law `(x1, y1) + (x2, y2) = (x3, y3)` with helper
/// columns `p = x1*x2`, `q = y1*y2`, `m = d*p*q`, transcribed from [`crate::add::EcAddAir::eval`].
///
/// `builder` is expected to be a selector-filtered sub-builder, which is exactly why this is a
/// free function over `AB: AirBuilder` rather than a `SubAirBuilder` embedding.
fn eval_add_law<AB: AirBuilder>(
    builder: &mut AB,
    d: AB::F,
    [x1, y1]: [AB::Var; 2],
    [x2, y2]: [AB::Var; 2],
    [x3, y3]: [AB::Var; 2],
    [p, q, m]: [AB::Var; 3],
) {
    // p = x1*x2, q = y1*y2
    builder.assert_zero(x1 * x2 - p);
    builder.assert_zero(y1 * y2 - q);

    // m = d*p*q
    builder.assert_zero((p * q) * d - m);

    // x3*(1+m) = x1*y2 + y1*x2
    builder.assert_zero(x3 * (AB::Expr::ONE + m) - (x1 * y2 + y1 * x2));

    // y3*(1-m) = q + p   (using a = -1, so y1*y2 - a*x1*x2 = q + p)
    builder.assert_zero(y3 * (AB::Expr::ONE - m) - (q + p));
}

impl<AB: AirBuilder> Air<AB> for EdDsaAir<AB::F>
where
    AB::F: Field,
{
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local_row = main.current_slice();
        let next_row = main.next_slice();

        // Selectors and the fixed-base doubling ladder. Copied out by value so the immutable
        // borrow of `builder` ends here.
        let &EdDsaPreprocessedCols {
            dbl_x,
            dbl_y,
            is_ladder,
            is_ladder_step,
            is_ladder_last,
            is_sgc_a,
            is_sgc_r,
            is_add,
        } = builder.preprocessed().current_slice().borrow();

        // Every row carries all four views of the same 45 columns; which one is meaningful is
        // decided by the selector each constraint below is gated on.
        let shared = *shared_cols(local_row);
        let next_shared = *shared_cols(next_row);
        let ladder = *ladder_cols(local_row);
        let next_ladder = *ladder_cols(next_row);
        let sgc = *sgc_cols(local_row);
        let add = *add_cols(local_row);
        let (var, fixed) = (ladder.var, ladder.fixed);
        let (next_var, next_fixed) = (next_ladder.var, next_ladder.fixed);

        // -- shared columns: constant down the whole table --------------------------------------
        //
        // This is the entire cross-row wiring mechanism: an AIR only sees `local`/`next`, so a
        // value produced on the ladder's last row and consumed three rows later has to ride along
        // on every row in between. Degree 1 each (degree 2 with the transition factor).
        {
            let mut transition = builder.when_transition();
            for (local_v, next_v) in [
                (shared.a_x, next_shared.a_x),
                (shared.a_y, next_shared.a_y),
                (shared.r_x, next_shared.r_x),
                (shared.r_y, next_shared.r_y),
                (shared.s, next_shared.s),
                (shared.k, next_shared.k),
                (shared.sb_x, next_shared.sb_x),
                (shared.sb_y, next_shared.sb_y),
                (shared.ka_x, next_shared.ka_x),
                (shared.ka_y, next_shared.ka_y),
            ] {
                transition.assert_eq(next_v, local_v);
            }
        }

        // -- public values, bound to the shared columns on the first row -------------------------
        //
        // Combined with the propagation above, this pins the shared columns on *every* row, which
        // is what makes them usable as the tie between row types.
        {
            let [a_x_pv, a_y_pv, r_x_pv, r_y_pv, s_pv, k_pv]: [AB::PublicVar;
                NUM_EDDSA_PUBLIC_VALUES] = builder.public_values().try_into().expect(
                "EdDsaAir expects exactly NUM_EDDSA_PUBLIC_VALUES public values: \
                 [a_x, a_y, r_x, r_y, s, k]",
            );
            let mut first_row = builder.when_first_row();
            first_row.assert_eq(shared.a_x, a_x_pv);
            first_row.assert_eq(shared.a_y, a_y_pv);
            first_row.assert_eq(shared.r_x, r_x_pv);
            first_row.assert_eq(shared.r_y, r_y_pv);
            first_row.assert_eq(shared.s, s_pv);
            first_row.assert_eq(shared.k, k_pv);
        }

        // -- ladder rows (`0 .. num_bits`): both double-and-add ladders in lockstep ---------------
        {
            let mut ladder_b = builder.when(is_ladder);

            // Variable-base ladder computing `[k]A`, transcribed from `EcScalarMulAir::eval`.
            // `acc + dbl` (the candidate next accumulator) and `dbl + dbl` (the doubling).
            eval_add_law(
                &mut ladder_b,
                self.d.dup(),
                [var.x_acc, var.y_acc],
                [var.x_dbl, var.y_dbl],
                [var.sum_x, var.sum_y],
                [var.add_p, var.add_q, var.add_m],
            );
            eval_add_law(
                &mut ladder_b,
                self.d.dup(),
                [var.dbl_x1, var.dbl_y1],
                [var.dbl_x2, var.dbl_y2],
                [var.dbl_x3, var.dbl_y3],
                [var.dbl_p, var.dbl_q, var.dbl_m],
            );
            // The doubling block doubles `(x_dbl, y_dbl)` with itself; wire both of its logical
            // inputs to the actual doubling-point columns.
            ladder_b.assert_eq(var.dbl_x1, var.x_dbl);
            ladder_b.assert_eq(var.dbl_y1, var.y_dbl);
            ladder_b.assert_eq(var.dbl_x2, var.x_dbl);
            ladder_b.assert_eq(var.dbl_y2, var.y_dbl);
            ladder_b.assert_bool(var.bit);
            ladder_b.assert_eq(
                var.result_x,
                var.bit * var.sum_x + (AB::Expr::ONE - var.bit) * var.x_acc,
            );
            ladder_b.assert_eq(
                var.result_y,
                var.bit * var.sum_y + (AB::Expr::ONE - var.bit) * var.y_acc,
            );

            // Fixed-base ladder computing `[s]B`, transcribed from `EcFixedScalarMulAir::eval`:
            // same shape, but `(x2, y2)` is the preprocessed `2^i * B` rather than a witness, so
            // there is no doubling block to constrain.
            eval_add_law(
                &mut ladder_b,
                self.d.dup(),
                [fixed.x_acc, fixed.y_acc],
                [dbl_x, dbl_y],
                [fixed.sum_x, fixed.sum_y],
                [fixed.add_p, fixed.add_q, fixed.add_m],
            );
            ladder_b.assert_bool(fixed.bit);
            ladder_b.assert_eq(
                fixed.result_x,
                fixed.bit * fixed.sum_x + (AB::Expr::ONE - fixed.bit) * fixed.x_acc,
            );
            ladder_b.assert_eq(
                fixed.result_y,
                fixed.bit * fixed.sum_y + (AB::Expr::ONE - fixed.bit) * fixed.y_acc,
            );
        }

        // -- ladder row 0: boundary values for both ladders --------------------------------------
        //
        // Row 0 is always a ladder row (`num_bits > 0`), so `is_first_row` alone is the right gate.
        {
            let mut first_row = builder.when_first_row();

            first_row.assert_zero(var.x_acc);
            first_row.assert_one(var.y_acc);
            first_row.assert_one(var.weight);
            first_row.assert_eq(var.scalar_acc, var.bit * var.weight);
            // Ties the variable ladder's base to the subgroup-checked `A`. `EcScalarMulAir`
            // leaves this to its caller — the base is a pure witness there — and it is what
            // closes the "scalar-mul base is never checked on-curve" gap.
            first_row.assert_eq(var.x_dbl, shared.a_x);
            first_row.assert_eq(var.y_dbl, shared.a_y);

            first_row.assert_zero(fixed.x_acc);
            first_row.assert_one(fixed.y_acc);
            first_row.assert_one(fixed.weight);
            first_row.assert_eq(fixed.scalar_acc, fixed.bit * fixed.weight);
            // The fixed ladder's base needs no constraint: `2^0 * B = B` holds by construction of
            // the preprocessed trace, which the prover cannot choose.
        }

        // -- ladder steps (`0 .. num_bits - 1`): advance both ladders ----------------------------
        //
        // A dedicated selector rather than `is_ladder * next.is_ladder`, so gating stays a single
        // degree-1 factor. It is zero on the table's last row, so no `when_transition` is needed.
        {
            let mut step = builder.when(is_ladder_step);

            step.assert_eq(
                next_var.x_acc,
                var.bit * var.sum_x + (AB::Expr::ONE - var.bit) * var.x_acc,
            );
            step.assert_eq(
                next_var.y_acc,
                var.bit * var.sum_y + (AB::Expr::ONE - var.bit) * var.y_acc,
            );
            step.assert_eq(next_var.x_dbl, var.dbl_x3);
            step.assert_eq(next_var.y_dbl, var.dbl_y3);
            step.assert_eq(next_var.weight, var.weight * AB::Expr::TWO);
            step.assert_eq(
                next_var.scalar_acc,
                var.scalar_acc + next_var.bit * next_var.weight,
            );

            step.assert_eq(
                next_fixed.x_acc,
                fixed.bit * fixed.sum_x + (AB::Expr::ONE - fixed.bit) * fixed.x_acc,
            );
            step.assert_eq(
                next_fixed.y_acc,
                fixed.bit * fixed.sum_y + (AB::Expr::ONE - fixed.bit) * fixed.y_acc,
            );
            step.assert_eq(next_fixed.weight, fixed.weight * AB::Expr::TWO);
            step.assert_eq(
                next_fixed.scalar_acc,
                fixed.scalar_acc + next_fixed.bit * next_fixed.weight,
            );
            // The fixed ladder's doubling point advances in the preprocessed trace, so there is
            // nothing to constrain here.
        }

        // -- ladder row `num_bits - 1`: read the results out into the shared columns --------------
        //
        // `scalar_acc` reaching `s`/`k` is what forces the ladders' bits to be the bit
        // decomposition of the public scalars (unique because `num_bits < F::bits()`).
        {
            let mut last = builder.when(is_ladder_last);
            last.assert_eq(fixed.scalar_acc, shared.s);
            last.assert_eq(fixed.result_x, shared.sb_x);
            last.assert_eq(fixed.result_y, shared.sb_y);
            last.assert_eq(var.scalar_acc, shared.k);
            last.assert_eq(var.result_x, shared.ka_x);
            last.assert_eq(var.result_y, shared.ka_y);
        }

        // -- subgroup-check rows (`num_bits`, `num_bits + 1`) ------------------------------------
        //
        // Transcribed from `EcSubgroupCheckAir::eval`: `(q_x, q_y) = 8P` for an on-curve witness
        // `P`, which certifies membership in the prime-order subgroup (cofactor 8). Identical on
        // both rows, so they share one selector; only the target differs, below.
        {
            let mut sgc_b = builder.when(is_sgc_a + is_sgc_r);

            // Feed P into the point-check block, then the curve equation
            // `-x^2 + y^2 = 1 + d*x^2*y^2` (transcribed from `EcPointCheckAir::eval`).
            sgc_b.assert_eq(sgc.cc_x, sgc.p_x);
            sgc_b.assert_eq(sgc.cc_y, sgc.p_y);
            sgc_b.assert_zero(sgc.cc_x * sgc.cc_x - sgc.cc_x2);
            sgc_b.assert_zero(sgc.cc_y * sgc.cc_y - sgc.cc_y2);
            sgc_b.assert_zero(
                (sgc.cc_y2 - sgc.cc_x2) - (AB::Expr::ONE + (sgc.cc_x2 * sgc.cc_y2) * self.d.dup()),
            );

            // Three chained doublings: P -> 2P -> 4P -> 8P.
            eval_add_law(
                &mut sgc_b,
                self.d.dup(),
                [sgc.d1_x1, sgc.d1_y1],
                [sgc.d1_x2, sgc.d1_y2],
                [sgc.d1_x3, sgc.d1_y3],
                [sgc.d1_p, sgc.d1_q, sgc.d1_m],
            );
            eval_add_law(
                &mut sgc_b,
                self.d.dup(),
                [sgc.d2_x1, sgc.d2_y1],
                [sgc.d2_x2, sgc.d2_y2],
                [sgc.d2_x3, sgc.d2_y3],
                [sgc.d2_p, sgc.d2_q, sgc.d2_m],
            );
            eval_add_law(
                &mut sgc_b,
                self.d.dup(),
                [sgc.d3_x1, sgc.d3_y1],
                [sgc.d3_x2, sgc.d3_y2],
                [sgc.d3_x3, sgc.d3_y3],
                [sgc.d3_p, sgc.d3_q, sgc.d3_m],
            );

            // Both logical inputs of each doubling block are the previous block's output.
            sgc_b.assert_eq(sgc.d1_x1, sgc.p_x);
            sgc_b.assert_eq(sgc.d1_y1, sgc.p_y);
            sgc_b.assert_eq(sgc.d1_x2, sgc.p_x);
            sgc_b.assert_eq(sgc.d1_y2, sgc.p_y);
            sgc_b.assert_eq(sgc.d2_x1, sgc.d1_x3);
            sgc_b.assert_eq(sgc.d2_y1, sgc.d1_y3);
            sgc_b.assert_eq(sgc.d2_x2, sgc.d1_x3);
            sgc_b.assert_eq(sgc.d2_y2, sgc.d1_y3);
            sgc_b.assert_eq(sgc.d3_x1, sgc.d2_x3);
            sgc_b.assert_eq(sgc.d3_y1, sgc.d2_y3);
            sgc_b.assert_eq(sgc.d3_x2, sgc.d2_x3);
            sgc_b.assert_eq(sgc.d3_y2, sgc.d2_y3);

            // The point under test is 8P.
            sgc_b.assert_eq(sgc.q_x, sgc.d3_x3);
            sgc_b.assert_eq(sgc.q_y, sgc.d3_y3);
        }

        // The target of each subgroup check, tied to the shared columns.
        {
            let mut sgc_a_b = builder.when(is_sgc_a);
            sgc_a_b.assert_eq(sgc.q_x, shared.a_x);
            sgc_a_b.assert_eq(sgc.q_y, shared.a_y);
        }
        {
            let mut sgc_r_b = builder.when(is_sgc_r);
            sgc_r_b.assert_eq(sgc.q_x, shared.r_x);
            sgc_r_b.assert_eq(sgc.q_y, shared.r_y);
        }

        // -- add row (`num_bits + 2`): `R + [k]A`, and the verification equation ------------------
        {
            let mut add_b = builder.when(is_add);
            eval_add_law(
                &mut add_b,
                self.d.dup(),
                [add.x1, add.y1],
                [add.x2, add.y2],
                [add.x3, add.y3],
                [add.p, add.q, add.m],
            );
            add_b.assert_eq(add.x1, shared.r_x);
            add_b.assert_eq(add.y1, shared.r_y);
            add_b.assert_eq(add.x2, shared.ka_x);
            add_b.assert_eq(add.y2, shared.ka_y);
            // `R + [k]A == [s]B`: the whole point of the circuit.
            add_b.assert_eq(add.x3, shared.sb_x);
            add_b.assert_eq(add.y3, shared.sb_y);
        }
    }
}
