//! QA encoder correctness: against the explicit generator matrix, against the
//! group-algebra definition, and the properties the PCS layers rely on.

use p3_brakedown::{LinearCode, LinearCodeExt};
use p3_field::PrimeCharacteristicRing;
use p3_goldilocks::Goldilocks;
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;
use p3_qa_code::{QaCode, generator_matrix, wht};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;

type F = Goldilocks;

fn code(log_n: usize, index: usize, seed: u64) -> QaCode<F> {
    QaCode::new_random(log_n, index, &mut ChaCha8Rng::seed_from_u64(seed))
}

fn rand_vec(rng: &mut ChaCha8Rng, n: usize) -> Vec<F> {
    (0..n).map(|_| rng.random::<F>()).collect()
}

/// `msg · G` with the explicit generator, for one message vector.
fn encode_by_generator(code: &QaCode<F>, msg: &[F]) -> Vec<F> {
    let g = generator_matrix(code);
    assert_eq!(g.height(), code.message_len());
    (0..g.width())
        .map(|col| {
            msg.iter()
                .enumerate()
                .map(|(r, &m)| m * g.values[r * g.width() + col])
                .sum()
        })
        .collect()
}

#[test]
fn dimensions() {
    for log_n in 1..6 {
        for index in [2, 3, 4] {
            let c = code(log_n, index, 7);
            assert_eq!(c.message_len(), 1 << log_n);
            assert_eq!(c.codeword_len(), index << log_n);
            assert_eq!(c.index(), index);
            assert_eq!(c.scaled_coeffs().len(), index - 1);
        }
    }
}

/// The fast encoder must equal multiplication by `(I | A_2 | … | A_c)`.
#[test]
fn matches_generator_matrix() {
    let mut rng = ChaCha8Rng::seed_from_u64(1);
    for log_n in 1..7 {
        for index in [2, 3, 4] {
            let c = code(log_n, index, 0xA1 + index as u64);
            let msg = rand_vec(&mut rng, c.message_len());
            assert_eq!(
                c.encode_vec(&msg),
                encode_by_generator(&c, &msg),
                "log_n={log_n} index={index}"
            );
        }
    }
}

/// Directly against the group-algebra definition: block `i` is the dyadic
/// convolution `(v · a_i)[k] = Σ_r v[r]·a_i[k ⊕ r]`.
#[test]
fn matches_dyadic_convolution() {
    let mut rng = ChaCha8Rng::seed_from_u64(2);
    let (log_n, index) = (5, 3);
    let c = code(log_n, index, 0xBEEF);
    let n = c.message_len();
    let msg = rand_vec(&mut rng, n);
    let codeword = c.encode_vec(&msg);

    for (i, scaled) in c.scaled_coeffs().iter().enumerate() {
        // a_i = WHT(2^{-n}·a_i').
        let mut a = scaled.clone();
        wht(&mut a);
        for k in 0..n {
            let expected: F = (0..n).map(|r| msg[r] * a[k ^ r]).sum();
            assert_eq!(codeword[(i + 1) * n + k], expected, "block {i} entry {k}");
        }
    }
}

/// Systematic: the message is the first block verbatim. Quasar leans on this
/// (`u₁ = v`, `q₁ = p`) to get its final evaluation claim for free.
#[test]
fn is_systematic() {
    let mut rng = ChaCha8Rng::seed_from_u64(3);
    let c = code(6, 4, 11);
    let msg = rand_vec(&mut rng, c.message_len());
    let codeword = c.encode_vec(&msg);
    assert_eq!(&codeword[..msg.len()], &msg[..]);
}

#[test]
fn is_linear() {
    let mut rng = ChaCha8Rng::seed_from_u64(4);
    let c = code(6, 2, 12);
    let n = c.message_len();
    let (x, y) = (rand_vec(&mut rng, n), rand_vec(&mut rng, n));
    let alpha: F = rng.random();

    let combo: Vec<F> = x.iter().zip(&y).map(|(&a, &b)| a * alpha + b).collect();
    let expected: Vec<F> = c
        .encode_vec(&x)
        .iter()
        .zip(c.encode_vec(&y).iter())
        .map(|(&a, &b)| a * alpha + b)
        .collect();
    assert_eq!(c.encode_vec(&combo), expected);
}

/// The batched matrix encoder (what the tensor PCS calls) must agree with the
/// per-column vector encoder.
#[test]
fn matrix_encode_matches_vector_encode() {
    let mut rng = ChaCha8Rng::seed_from_u64(5);
    for (log_n, index, width) in [(4, 2, 1), (5, 2, 3), (6, 3, 8), (7, 4, 5)] {
        let c = code(log_n, index, 0x1234);
        let (n, cw) = (c.message_len(), c.codeword_len());
        let values = rand_vec(&mut rng, n * width);
        let msg = RowMajorMatrix::new(values.clone(), width);

        let got = c.encode(&msg);
        assert_eq!(got.height(), cw);
        assert_eq!(got.width(), width);

        for col in 0..width {
            let column: Vec<F> = (0..n).map(|r| values[r * width + col]).collect();
            let expected = c.encode_vec(&column);
            for (r, e) in expected.iter().enumerate() {
                assert_eq!(got.values[r * width + col], *e, "w={width} col={col} row={r}");
            }
        }
    }
}

/// `encode_ext` must apply the same `F`-linear code to extension-field
/// messages: the tensor PCS commits base-field rows but encodes the combined
/// row over the challenge field.
#[test]
fn encode_ext_agrees_on_embedded_values() {
    use p3_goldilocks::Goldilocks as Base;
    type EF = p3_field::extension::BinomialExtensionField<Base, 2>;

    let mut rng = ChaCha8Rng::seed_from_u64(6);
    let c = code(5, 2, 0x5EED);
    let n = c.message_len();
    let width = 3;
    let base = rand_vec(&mut rng, n * width);

    let base_encoded = c.encode(&RowMajorMatrix::new(base.clone(), width));
    let lifted: Vec<EF> = base.iter().map(|&x| EF::from(x)).collect();
    let ext_encoded = <QaCode<F> as LinearCodeExt<F, EF>>::encode_ext(
        &c,
        &RowMajorMatrix::new(lifted, width),
    );

    for (b, e) in base_encoded.values.iter().zip(ext_encoded.values.iter()) {
        assert_eq!(EF::from(*b), *e);
    }
}

/// `middle_vec` is the encoding sumcheck's witness `v' = WHT(v)`; the parity
/// blocks must be recoverable from it, since that is exactly what the
/// code-switching argument certifies.
#[test]
fn middle_vec_reproduces_parity_blocks() {
    let mut rng = ChaCha8Rng::seed_from_u64(7);
    let c = code(6, 3, 0xC0DE);
    let n = c.message_len();
    let msg = rand_vec(&mut rng, n);

    let middle = c.middle_vec(&msg);
    let codeword = c.encode_vec(&msg);

    for (i, coeff) in c.scaled_coeffs().iter().enumerate() {
        let mut block: Vec<F> = middle.iter().zip(coeff).map(|(&m, &e)| m * e).collect();
        wht(&mut block);
        assert_eq!(&codeword[(i + 1) * n..(i + 2) * n], &block[..], "block {i}");
    }
}

/// A weight-1 message difference is the worst case for a *systematic* code: the
/// message block can only differ in one position, so all the distance has to
/// come from the parity blocks. Sanity-checks the code is not degenerate (the
/// real distance bound lives in `p3-security`); a zero or sparse coefficient
/// vector would fail this.
#[test]
fn weight_one_difference_spreads_across_the_parity_blocks() {
    let mut rng = ChaCha8Rng::seed_from_u64(8);
    let c = code(6, 2, 0xD15);
    let n = c.message_len();
    let x = rand_vec(&mut rng, n);
    let mut y = x.clone();
    y[3] += F::ONE;

    let (cx, cy) = (c.encode_vec(&x), c.encode_vec(&y));

    // Every parity position must move.
    let parity_agree = cx[n..]
        .iter()
        .zip(&cy[n..])
        .filter(|(a, b)| a == b)
        .count();
    assert_eq!(parity_agree, 0, "parity block agreed in {parity_agree}/{n}");

    // Hence relative distance ≈ 1/2 — the ceiling for a systematic rate-1/2
    // code at message weight 1, and comfortably above the proven δ ≈ 0.46.
    let differ = cx.iter().zip(&cy).filter(|(a, b)| a != b).count();
    assert!(
        differ * 100 >= 45 * cx.len(),
        "relative distance {differ}/{} too small",
        cx.len()
    );
}
