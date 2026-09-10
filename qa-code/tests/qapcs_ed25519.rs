//! **QAPCS**: the Brakedown-style tensor PCS instantiated with a QA code, over
//! the native ed25519 base field (`p = 2²⁵⁵ − 19`, non-two-adic).
//!
//! This is the baseline Quasar is measured against, and it comes for free:
//! `QaCode` implements `LinearCodeExt`, so `TensorPcs` accepts it unchanged.
//! The query count comes from the proven QA distance in `p3-security`, not from
//! a hand-picked constant.

use p3_brakedown::{LinearCode, TensorPcs, TensorPcsError};
use p3_challenger::DuplexChallengerPF;
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519, Poseidon2Ed25519,
};
use p3_field::PrimeCharacteristicRing;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_multilinear_util::point::Point;
use p3_qa_code::QaCode;
use p3_security::QaRegime;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;

type F = Ed25519;
type EF = Ed25519;

const WIDTH: usize = 3;
const RATE: usize = 2;
const OUT: usize = 1;

type Perm = Poseidon2Ed25519<WIDTH>;
type MyHash = PaddingFreeSponge<Perm, WIDTH, RATE, OUT>;
type MyCompress = TruncatedPermutation<Perm, 2, OUT, WIDTH>;
type MyMmcs = MerkleTreeMmcs<Ed25519, Ed25519, MyHash, MyCompress, 2, OUT>;
type Challenger = DuplexChallengerPF<Ed25519, Perm, WIDTH, RATE>;
type Pcs = TensorPcs<F, EF, QaCode<F>, MyMmcs>;

const LOG_ROWS: usize = 3; // 8 rows
const LOG_COLS: usize = 6; // QA message length 64
const NUM_VARS: usize = LOG_ROWS + LOG_COLS;
const INDEX: usize = 2; // rate 1/2
/// The instance is far too small for the real query count (~420 columns > the
/// 128-column codeword), so tests use a token count; the *real* count is
/// asserted separately in [`real_query_count_comes_from_the_distance_bound`].
const NUM_QUERIES: usize = 16;

fn perm() -> Perm {
    let mut rng = ChaCha8Rng::seed_from_u64(0xC0FFEE);
    Poseidon2Ed25519::<WIDTH>::new_from_rng(
        ED25519_POSEIDON2_FULL_ROUNDS,
        ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
        &mut rng,
    )
}

fn challenger() -> Challenger {
    Challenger::new(perm())
}

fn pcs() -> Pcs {
    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);
    let mut code_rng = ChaCha8Rng::seed_from_u64(0x0A_C0DE);
    let code = QaCode::<F>::new_random(LOG_COLS, INDEX, &mut code_rng);
    assert_eq!(code.message_len(), 1 << LOG_COLS);
    assert_eq!(code.codeword_len(), INDEX << LOG_COLS);
    assert!(NUM_QUERIES <= code.codeword_len());
    TensorPcs::new(code, mmcs, LOG_ROWS, NUM_QUERIES)
}

fn rand_evals(rng: &mut ChaCha8Rng) -> Vec<F> {
    (0..(1 << NUM_VARS)).map(|_| rng.random::<F>()).collect()
}

fn rand_point(rng: &mut ChaCha8Rng) -> Point<EF> {
    Point::new((0..NUM_VARS).map(|_| rng.random::<EF>()).collect())
}

/// Independent multilinear evaluation, big-endian (coordinate 0 = most
/// significant bit) — the PCS convention.
fn multilinear_eval(evals: &[F], z: &Point<EF>) -> EF {
    let mut table = vec![EF::ONE];
    for &zi in z.as_slice() {
        let mut next = Vec::with_capacity(table.len() * 2);
        for &t in &table {
            next.push(t * (EF::ONE - zi));
            next.push(t * zi);
        }
        table = next;
    }
    evals.iter().zip(&table).map(|(&e, &eq)| eq * e).sum()
}

#[test]
fn qapcs_roundtrip_and_matches_direct_eval() {
    let mut rng = ChaCha8Rng::seed_from_u64(1);
    let pcs = pcs();
    let evals = rand_evals(&mut rng);
    let point = rand_point(&mut rng);

    let mut ch_p = challenger();
    let (commitment, prover_data) = pcs.commit(evals.clone(), &mut ch_p);
    let (value, proof) = pcs.open(&prover_data, &point, &mut ch_p);

    assert_eq!(value, multilinear_eval(&evals, &point), "opened value wrong");

    let mut ch_v = challenger();
    pcs.verify(&commitment, &point, value, &proof, &mut ch_v)
        .expect("honest QAPCS opening must verify over ed25519");
}

#[test]
fn qapcs_rejects_tampered_column() {
    let mut rng = ChaCha8Rng::seed_from_u64(3);
    let pcs = pcs();
    let evals = rand_evals(&mut rng);
    let point = rand_point(&mut rng);

    let mut ch_p = challenger();
    let (commitment, prover_data) = pcs.commit(evals, &mut ch_p);
    let (value, mut proof) = pcs.open(&prover_data, &point, &mut ch_p);

    proof.openings[0][0][0] += F::ONE;

    let mut ch_v = challenger();
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut ch_v)
        .expect_err("a tampered column opening must be rejected");
    assert!(matches!(err, TensorPcsError::Merkle(_)), "got {err:?}");
}

#[test]
fn qapcs_rejects_wrong_value() {
    let mut rng = ChaCha8Rng::seed_from_u64(2);
    let pcs = pcs();
    let evals = rand_evals(&mut rng);
    let point = rand_point(&mut rng);

    let mut ch_p = challenger();
    let (commitment, prover_data) = pcs.commit(evals, &mut ch_p);
    let (value, proof) = pcs.open(&prover_data, &point, &mut ch_p);

    let mut ch_v = challenger();
    let err = pcs
        .verify(&commitment, &point, value + EF::ONE, &proof, &mut ch_v)
        .expect_err("a wrong claimed value must be rejected");
    assert!(
        matches!(err, TensorPcsError::ClaimedValueMismatch),
        "got {err:?}"
    );
}

/// Rate 1/4 works the same way and needs fewer columns for the same security.
#[test]
fn qapcs_roundtrip_at_index_four() {
    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);
    let mut code_rng = ChaCha8Rng::seed_from_u64(0x4444);
    let code = QaCode::<F>::new_random(LOG_COLS, 4, &mut code_rng);
    let pcs: Pcs = TensorPcs::new(code, mmcs, LOG_ROWS, NUM_QUERIES);

    let mut rng = ChaCha8Rng::seed_from_u64(9);
    let evals = rand_evals(&mut rng);
    let point = rand_point(&mut rng);

    let mut ch_p = challenger();
    let (commitment, prover_data) = pcs.commit(evals.clone(), &mut ch_p);
    let (value, proof) = pcs.open(&prover_data, &point, &mut ch_p);
    assert_eq!(value, multilinear_eval(&evals, &point));

    let mut ch_v = challenger();
    pcs.verify(&commitment, &point, value, &proof, &mut ch_v)
        .expect("index-4 QAPCS opening must verify");
}

/// The number of columns a *real* QAPCS instance must open, from the proven QA
/// distance over ed25519. Documents why the tests above use a token count.
#[test]
fn real_query_count_comes_from_the_distance_bound() {
    let regime = QaRegime {
        row_log_size: 16,
        index: INDEX,
        field_bits: 255,
        distance_failure_bits: 103,
        query_security_bits: 103,
    };
    assert_eq!(regime.num_queries(), Some(425));
    // Far more than this test's 128-element codeword, hence the token count.
    assert!(regime.num_queries().unwrap() > (INDEX << LOG_COLS));
}
