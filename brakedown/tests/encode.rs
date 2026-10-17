//! Encoder correctness tests for the Brakedown linear code.
//!
//! These validate the *mechanics* (M2 scope): every code here is a genuine
//! linear code, the Brakedown construction is systematic and has the right
//! dimensions, and the RS base has RS distance. Production parameters and
//! soundness are out of scope until M4.

use p3_baby_bear::BabyBear;
use p3_brakedown::{BrakedownCode, LinearCode, LinearCodeExt, ReedSolomonCode, SparseMatrix};
use p3_field::PrimeCharacteristicRing;
use p3_field::extension::BinomialExtensionField;
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;
use rand::SeedableRng;
use rand::rngs::SmallRng;

type F = BabyBear;
type EF = BinomialExtensionField<F, 4>;

/// Deterministic RNG so every run exercises the same random code + messages.
fn rng() -> SmallRng {
    SmallRng::seed_from_u64(0xB2A4E_D0)
}

/// A random `rows × width` message matrix.
fn rand_message(rng: &mut SmallRng, rows: usize, width: usize) -> RowMajorMatrix<F> {
    use rand::RngExt;
    let values = (0..rows * width).map(|_| rng.random::<F>()).collect();
    RowMajorMatrix::new(values, width)
}

/// `alpha * a + beta * b`, elementwise (a, b same shape).
fn lin_combine(a: &RowMajorMatrix<F>, b: &RowMajorMatrix<F>, alpha: F, beta: F) -> RowMajorMatrix<F> {
    assert_eq!(a.values.len(), b.values.len());
    let values = a
        .values
        .iter()
        .zip(&b.values)
        .map(|(&x, &y)| alpha * x + beta * y)
        .collect();
    RowMajorMatrix::new(values, a.width())
}

/// Assert `code` is linear: encode(αx + βy) == α·encode(x) + β·encode(y).
fn assert_linear<C: LinearCode<F>>(rng: &mut SmallRng, code: &C, width: usize) {
    use rand::RngExt;
    let x = rand_message(rng, code.message_len(), width);
    let y = rand_message(rng, code.message_len(), width);
    let alpha: F = rng.random();
    let beta: F = rng.random();

    let lhs = code.encode(&lin_combine(&x, &y, alpha, beta));
    let rhs = lin_combine(&code.encode(&x), &code.encode(&y), alpha, beta);
    assert_eq!(lhs.values, rhs.values, "code is not linear");
}

#[test]
fn reed_solomon_dims_linearity_and_distance() {
    let mut rng = rng();
    let (k, n) = (6, 12);
    let rs = ReedSolomonCode::<F>::new(k, n);
    assert_eq!(rs.message_len(), k);
    assert_eq!(rs.codeword_len(), n);

    // Dimensions of an actual encoding.
    let msg = rand_message(&mut rng, k, 3);
    let cw = rs.encode(&msg);
    assert_eq!(cw.height(), n);
    assert_eq!(cw.width(), 3);

    assert_linear(&mut rng, &rs, 4);

    // RS distance: two messages differing in one coordinate differ in >= n-k+1
    // codeword positions. Use single-column messages so "position" = row.
    let m0 = rand_message(&mut rng, k, 1);
    let mut m1 = m0.clone();
    m1.values[0] += F::ONE; // perturb one message coordinate
    let c0 = rs.encode(&m0);
    let c1 = rs.encode(&m1);
    let differing = c0
        .values
        .iter()
        .zip(&c1.values)
        .filter(|(a, b)| a != b)
        .count();
    assert!(
        differing >= n - k + 1,
        "RS distance violated: {differing} < {}",
        n - k + 1
    );
}

/// Build a one-level Brakedown code: message_len 16, RS base (k=8, n=12).
///
///   A: 8 × 16  (message → inner input)
///   inner RS: k=8, n=12
///   B: 6 × 12  (inner codeword → tail)
///   codeword_len = 16 + 12 + 6 = 34
fn one_level(rng: &mut SmallRng) -> BrakedownCode<F, EF> {
    let a = SparseMatrix::<F>::rand_fixed_row_weight(rng, 8, 16, 4);
    let inner = ReedSolomonCode::<F>::new(8, 12);
    let b = SparseMatrix::<F>::rand_fixed_row_weight(rng, 6, 12, 3);
    BrakedownCode::new(a, b, Box::new(inner))
}

#[test]
fn brakedown_one_level_dims_systematic_and_linear() {
    let mut rng = rng();
    let code = one_level(&mut rng);
    assert_eq!(code.message_len(), 16);
    assert_eq!(code.codeword_len(), 34);

    let width = 3;
    let msg = rand_message(&mut rng, 16, width);
    let cw = code.encode(&msg);
    assert_eq!(cw.height(), 34);
    assert_eq!(cw.width(), width);

    // Systematic: the codeword's top `message_len` rows equal the message.
    assert_eq!(
        &cw.values[..16 * width],
        &msg.values[..],
        "Brakedown code is not systematic"
    );

    assert_linear(&mut rng, &code, width);
}

#[test]
fn brakedown_two_level_dims_systematic_and_linear() {
    // Level 1 (inner): message_len 16 → codeword_len 34 (as above).
    // Level 0 (outer): message_len 32.
    //   A: 16 × 32 ; inner = level-1 Brakedown (message_len 16, codeword_len 34)
    //   B: 12 × 34 ; codeword_len = 32 + 34 + 12 = 78
    let mut rng = rng();
    let inner = one_level(&mut rng);
    assert_eq!(inner.message_len(), 16);
    assert_eq!(inner.codeword_len(), 34);

    let a = SparseMatrix::<F>::rand_fixed_row_weight(&mut rng, 16, 32, 5);
    let b = SparseMatrix::<F>::rand_fixed_row_weight(&mut rng, 12, 34, 4);
    let code = BrakedownCode::new(a, b, Box::new(inner));

    assert_eq!(code.message_len(), 32);
    assert_eq!(code.codeword_len(), 78);

    let width = 2;
    let msg = rand_message(&mut rng, 32, width);
    let cw = code.encode(&msg);
    assert_eq!(cw.height(), 78);
    assert_eq!(&cw.values[..32 * width], &msg.values[..], "outer code not systematic");

    assert_linear(&mut rng, &code, width);
}

/// A random `rows × width` message matrix over `EF`.
fn rand_message_ef(rng: &mut SmallRng, rows: usize, width: usize) -> RowMajorMatrix<EF> {
    use rand::RngExt;
    let values = (0..rows * width).map(|_| rng.random::<EF>()).collect();
    RowMajorMatrix::new(values, width)
}

/// `alpha * a + beta * b`, elementwise, over `EF`.
fn lin_combine_ef(
    a: &RowMajorMatrix<EF>,
    b: &RowMajorMatrix<EF>,
    alpha: EF,
    beta: EF,
) -> RowMajorMatrix<EF> {
    let values = a
        .values
        .iter()
        .zip(&b.values)
        .map(|(&x, &y)| alpha * x + beta * y)
        .collect();
    RowMajorMatrix::new(values, a.width())
}

#[test]
fn brakedown_encode_ext_is_linear_and_agrees_with_base() {
    use rand::RngExt;
    let mut rng = rng();
    let code = one_level(&mut rng);
    let k = code.message_len();
    let width = 3;

    // `encode_ext` is EF-linear.
    let x = rand_message_ef(&mut rng, k, width);
    let y = rand_message_ef(&mut rng, k, width);
    let alpha: EF = rng.random();
    let beta: EF = rng.random();
    let lhs = code.encode_ext(&lin_combine_ef(&x, &y, alpha, beta));
    let rhs = lin_combine_ef(&code.encode_ext(&x), &code.encode_ext(&y), alpha, beta);
    assert_eq!(lhs.values, rhs.values, "encode_ext is not linear");

    // When the EF message is base-field-valued, encode_ext agrees with encode
    // (same F-linear generator, lifted componentwise).
    let base = rand_message(&mut rng, k, width);
    let base_as_ef = RowMajorMatrix::new(
        base.values.iter().map(|&v| EF::from(v)).collect(),
        width,
    );
    let cw_base = code.encode(&base);
    let cw_ext = code.encode_ext(&base_as_ef);
    let cw_base_lifted: Vec<EF> = cw_base.values.iter().map(|&v| EF::from(v)).collect();
    assert_eq!(
        cw_ext.values, cw_base_lifted,
        "encode_ext must agree with encode on base-field-valued messages"
    );
}

#[test]
fn brakedown_from_message_len_recurses_systematic_linear_and_ext() {
    // A power-of-two top message length that recurses at least one level
    // (32 > base_threshold 16 → inner m = round(0.238·32) = 8 → RS base).
    let mut rng = rng();
    let message_len = 32;
    let code = BrakedownCode::<F, EF>::from_message_len(message_len, &mut rng);
    assert_eq!(code.message_len(), message_len);
    let n = code.codeword_len();
    // Rate ≈ 1.72 ⇒ codeword clearly longer than the message, but not doubled.
    assert!(n > message_len, "codeword must be longer than the message");
    assert!(n < 3 * message_len, "codeword length {n} unexpectedly large");

    let width = 3;
    let msg = rand_message(&mut rng, message_len, width);
    let cw = code.encode(&msg);
    assert_eq!(cw.height(), n);
    assert_eq!(cw.width(), width);
    // Systematic.
    assert_eq!(
        &cw.values[..message_len * width],
        &msg.values[..],
        "from_message_len code is not systematic"
    );

    assert_linear(&mut rng, &code, width);

    // encode_ext linear + agrees with encode on base-valued messages.
    use rand::RngExt;
    let x = rand_message_ef(&mut rng, message_len, width);
    let y = rand_message_ef(&mut rng, message_len, width);
    let alpha: EF = rng.random();
    let beta: EF = rng.random();
    let lhs = code.encode_ext(&lin_combine_ef(&x, &y, alpha, beta));
    let rhs = lin_combine_ef(&code.encode_ext(&x), &code.encode_ext(&y), alpha, beta);
    assert_eq!(lhs.values, rhs.values, "recursive encode_ext is not linear");

    let base_as_ef =
        RowMajorMatrix::new(msg.values.iter().map(|&v| EF::from(v)).collect(), width);
    let cw_ext = code.encode_ext(&base_as_ef);
    let cw_base_lifted: Vec<EF> = cw.values.iter().map(|&v| EF::from(v)).collect();
    assert_eq!(
        cw_ext.values, cw_base_lifted,
        "recursive encode_ext must agree with encode on base-valued messages"
    );
}

#[test]
#[should_panic]
fn brakedown_new_rejects_mismatched_inner_input() {
    // A.height() (7) must equal inner.message_len() (8) — this must panic.
    let mut rng = rng();
    let a = SparseMatrix::<F>::rand_fixed_row_weight(&mut rng, 7, 16, 4);
    let inner = ReedSolomonCode::<F>::new(8, 12);
    let b = SparseMatrix::<F>::rand_fixed_row_weight(&mut rng, 6, 12, 3);
    let _: BrakedownCode<F, EF> = BrakedownCode::new(a, b, Box::new(inner));
}
