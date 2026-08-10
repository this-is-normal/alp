use num_bigint::BigUint;
use p3_bn254::Bn254;
use p3_ed25519_air::add::add_points;
use p3_ed25519_air::eddsa::{
    NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS, add_row, generate_preprocessed_trace,
    generate_trace_rows, ladder_cols, ladder_results, public_values, sgc_a_row, sgc_r_row,
    shared_cols, trace_height,
};
use p3_field::PrimeCharacteristicRing;
use p3_matrix::Matrix;

#[test]
fn smoke() {
    type F = Bn254;
    let d = F::from_u64(7);
    let identity = (F::ZERO, F::ONE);
    let base_b = identity;
    let l = BigUint::from(11u32);
    let n = 5;
    let s_bits = vec![false; n];
    let k_bits = vec![true, false, true, true, false];

    let a = identity;
    let r = identity;

    let trace = generate_trace_rows(d, base_b, &l, a, r, &s_bits, &k_bits);
    let pre = generate_preprocessed_trace(d, base_b, n);
    assert_eq!(trace.height(), trace_height(n));
    assert_eq!(pre.height(), trace_height(n));
    assert_eq!(trace.width(), NUM_EDDSA_COLS);
    assert_eq!(pre.width(), NUM_EDDSA_PREPROCESSED_COLS);

    let pv = public_values(a, r, &s_bits, &k_bits);
    // row-0 shared columns match the public values
    let row0 = shared_cols(&trace.values[..NUM_EDDSA_COLS]);
    assert_eq!(
        vec![row0.a_x, row0.a_y, row0.r_x, row0.r_y, row0.s, row0.k],
        pv
    );

    // shared columns constant on every row
    for row in 0..trace.height() {
        let c = shared_cols(&trace.values[row * NUM_EDDSA_COLS..(row + 1) * NUM_EDDSA_COLS]);
        assert_eq!((c.a_x, c.a_y), a);
        assert_eq!(c.k, row0.k);
        assert_eq!((c.ka_x, c.ka_y), (row0.ka_x, row0.ka_y));
    }

    // ladder last row results match the shared sb / ka
    let (sb, ka) = ladder_results(&trace, n);
    assert_eq!(sb, (row0.sb_x, row0.sb_y));
    assert_eq!(ka, (row0.ka_x, row0.ka_y));

    // ladder row 0: acc = identity, weight = 1, var base = A
    let l0 = ladder_cols(&trace.values[..NUM_EDDSA_COLS]);
    assert_eq!((l0.var.x_acc, l0.var.y_acc), identity);
    assert_eq!((l0.var.x_dbl, l0.var.y_dbl), a);
    assert_eq!(l0.var.weight, F::ONE);
    assert_eq!(l0.fixed.weight, F::ONE);
    assert_eq!((l0.fixed.x_acc, l0.fixed.y_acc), identity);

    // preprocessed dbl ladder feeds the fixed block: dbl_i = 2^i * B
    let mut cur = base_b;
    for row in 0..n {
        let p = &pre.values[row * NUM_EDDSA_PREPROCESSED_COLS..];
        assert_eq!((p[0], p[1]), cur);
        assert_eq!(p[2], F::ONE, "is_ladder");
        cur = add_points(d, cur, cur);
    }
    // selectors on the aux rows
    let sel = |row: usize, i: usize| pre.values[row * NUM_EDDSA_PREPROCESSED_COLS + i];
    assert_eq!(sel(n - 1, 4), F::ONE, "is_ladder_last");
    assert_eq!(sel(sgc_a_row(n), 5), F::ONE);
    assert_eq!(sel(sgc_r_row(n), 6), F::ONE);
    assert_eq!(sel(add_row(n), 7), F::ONE);
    for row in add_row(n) + 1..trace.height() {
        for i in 2..NUM_EDDSA_PREPROCESSED_COLS {
            assert_eq!(sel(row, i), F::ZERO, "padding row {row} selector {i}");
        }
    }
}
