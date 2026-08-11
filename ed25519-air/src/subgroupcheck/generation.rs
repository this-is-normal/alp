use core::borrow::BorrowMut;

use alloc::vec::Vec;

use num_bigint::BigUint;
use p3_field::Field;
use p3_matrix::dense::RowMajorMatrix;

use crate::add::{EcAddCols, Point, add_points};
use crate::subgroupcheck::columns::{EcSubgroupCheckCols, NUM_EC_SUBGROUPCHECK_COLS};

/// Doubles `pt` (`a = -1`), returning the filled `EcAddCols` for that step (`x1,y1,x2,y2` are
/// both copies of `pt`; `x3,y3` is the result).
fn double_step<F: Field>(d: F, pt: Point<F>) -> EcAddCols<F> {
    let (x, y) = pt;
    let p = x * x;
    let q = y * y;
    let m = d * p * q;
    let (x3, y3) = add_points(d, pt, pt);
    EcAddCols {
        x1: x,
        y1: y,
        x2: x,
        y2: y,
        x3,
        y3,
        p,
        q,
        m,
    }
}

/// Generates a trace with one row per witness point `P`, filling in `8*P` and every helper
/// witness column along the way.
///
/// Pads with the identity `(0, 1)`, which is on-curve for any `a, d` and its own double, so
/// all three doubling blocks stay valid on padding rows.
pub fn generate_trace_rows<F: Field>(d: F, witnesses: &[(F, F)]) -> RowMajorMatrix<F> {
    let num_rows = witnesses.len().next_power_of_two().max(1);
    let mut trace = RowMajorMatrix::new(
        F::zero_vec(num_rows * NUM_EC_SUBGROUPCHECK_COLS),
        NUM_EC_SUBGROUPCHECK_COLS,
    );

    let identity = (F::ZERO, F::ONE);
    for (row, chunk) in trace
        .values
        .chunks_mut(NUM_EC_SUBGROUPCHECK_COLS)
        .enumerate()
    {
        let cols: &mut EcSubgroupCheckCols<F> = chunk.borrow_mut();
        let (p_x, p_y) = witnesses.get(row).copied().unwrap_or(identity);

        let cc_x2 = p_x * p_x;
        let cc_y2 = p_y * p_y;

        let d1 = double_step(d, (p_x, p_y));
        let d2 = double_step(d, (d1.x3, d1.y3));
        let d3 = double_step(d, (d2.x3, d2.y3));

        *cols = EcSubgroupCheckCols {
            p_x,
            p_y,
            cc_x: p_x,
            cc_y: p_y,
            cc_x2,
            cc_y2,
            d1_x1: d1.x1,
            d1_y1: d1.y1,
            d1_x2: d1.x2,
            d1_y2: d1.y2,
            d1_x3: d1.x3,
            d1_y3: d1.y3,
            d1_p: d1.p,
            d1_q: d1.q,
            d1_m: d1.m,
            d2_x1: d2.x1,
            d2_y1: d2.y1,
            d2_x2: d2.x2,
            d2_y2: d2.y2,
            d2_x3: d2.x3,
            d2_y3: d2.y3,
            d2_p: d2.p,
            d2_q: d2.q,
            d2_m: d2.m,
            d3_x1: d3.x1,
            d3_y1: d3.y1,
            d3_x2: d3.x2,
            d3_y2: d3.y2,
            d3_x3: d3.x3,
            d3_y3: d3.y3,
            d3_p: d3.p,
            d3_q: d3.q,
            d3_m: d3.m,
            q_x: d3.x3,
            q_y: d3.y3,
        };
    }

    trace
}

/// Multiplies `pt` by `scalar` via double-and-add, reusing `EcAddAir`'s complete addition law
/// (so it's valid at every step, including when the accumulator or base point is the identity).
pub fn scalar_mul<F: Field>(d: F, scalar: &BigUint, pt: Point<F>) -> Point<F> {
    let mut acc = (F::ZERO, F::ONE);
    let mut base = pt;
    for i in 0..scalar.bits() {
        if scalar.bit(i) {
            acc = add_points(d, acc, base);
        }
        base = add_points(d, base, base);
    }
    acc
}

/// Modular inverse of `8` mod prime `l`, via Fermat's little theorem (`8^(l-2) mod l`).
///
/// Requires `l` to be prime and coprime to 8 (true of ed25519's subgroup order), since
/// Fermat's little theorem is what makes `8 * inverse_of_8_mod(l) == 1 (mod l)` hold.
fn inverse_of_8_mod(l: &BigUint) -> BigUint {
    BigUint::from(8u32).modpow(&(l - BigUint::from(2u32)), l)
}

/// Derives the witness `P = (8^-1 mod l) * Q` for a target point `Q`, given the known order
/// `l` of the prime-order subgroup.
///
/// `8*P == Q` iff `Q`'s order divides `l` (Lagrange) — i.e. iff `Q` genuinely belongs to that
/// subgroup. For an in-subgroup `Q` this always produces a valid membership witness; for a
/// `Q` outside the subgroup it produces *some* point, but `8*P != Q`, so
/// `EcSubgroupCheckAir`'s boundary constraint (still) correctly rejects it.
pub fn witness_for_target<F: Field>(d: F, l: &BigUint, target: Point<F>) -> Point<F> {
    scalar_mul(d, &inverse_of_8_mod(l), target)
}

/// Like [`generate_trace_rows`], but takes the target points `Q` directly and derives each
/// row's witness `P` internally via [`witness_for_target`], given the subgroup's known order
/// `l`. This removes the need for callers to solve for `P` themselves outside the circuit.
pub fn generate_trace_rows_for_targets<F: Field>(
    d: F,
    l: &BigUint,
    targets: &[(F, F)],
) -> RowMajorMatrix<F> {
    let witnesses: Vec<Point<F>> = targets
        .iter()
        .map(|&target| witness_for_target(d, l, target))
        .collect();
    generate_trace_rows(d, &witnesses)
}
