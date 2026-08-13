//! Tensor-PCS tests using the **real recursive Brakedown code over the native
//! ed25519 base field** (`p = 2^255 − 19`, non-two-adic, non-`PrimeField64`).
//!
//! This is the ed25519 counterpart to `tensor_pcs_brakedown.rs` (BabyBear). It
//! validates the Brakedown tensor commitment *directly* on the target field —
//! commit/open/verify, the opened value matching an independent multilinear
//! evaluation, and rejection of a tampered column — with the width-3
//! `Poseidon2Ed25519` Merkle stack and the general-`PrimeField` duplex challenger.

use p3_brakedown::{BrakedownCode, LinearCode, TensorPcs, TensorPcsError};
use p3_challenger::DuplexChallengerPF;
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519, Poseidon2Ed25519,
};
use p3_field::PrimeCharacteristicRing;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_multilinear_util::point::Point;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

type F = Ed25519;
type EF = Ed25519;

// ed25519 Poseidon2 supports width 3; capacity 1, digest 1 element (~127-bit collisions).
const WIDTH: usize = 3;
const RATE: usize = 2;
const OUT: usize = 1;

type Perm = Poseidon2Ed25519<WIDTH>;
type MyHash = PaddingFreeSponge<Perm, WIDTH, RATE, OUT>;
type MyCompress = TruncatedPermutation<Perm, 2, OUT, WIDTH>;
type MyMmcs = MerkleTreeMmcs<Ed25519, Ed25519, MyHash, MyCompress, 2, OUT>;
type Challenger = DuplexChallengerPF<Ed25519, Perm, WIDTH, RATE>;
type Pcs = TensorPcs<F, EF, BrakedownCode<F, EF>, MyMmcs>;

const LOG_ROWS: usize = 3; // rows = 8
const COLS: usize = 32; // Brakedown message length (power of two; recurses one level)
const NUM_QUERIES: usize = 12;
const NUM_VARS: usize = 8; // rows(8) * cols(32) = 256 = 2^8

fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(0xC0FFEE);
    Poseidon2Ed25519::<WIDTH>::new_from_rng(
        ED25519_POSEIDON2_FULL_ROUNDS,
        ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
        &mut rng,
    )
}

fn challenger() -> Challenger {
    Challenger::new(perm())
}

/// Build the PCS with a real Brakedown code, drawn from a fixed-seed RNG so both
/// prover and verifier share one code.
fn pcs() -> Pcs {
    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);
    let mut code_rng = SmallRng::seed_from_u64(0xBADC0DE);
    let code = BrakedownCode::<F, EF>::from_message_len(COLS, &mut code_rng);
    assert_eq!(code.message_len(), COLS);
    assert!(NUM_QUERIES <= code.codeword_len());
    TensorPcs::new(code, mmcs, LOG_ROWS, NUM_QUERIES)
}

fn rand_evals(rng: &mut SmallRng) -> Vec<F> {
    (0..(1 << NUM_VARS)).map(|_| rng.random::<F>()).collect()
}

fn rand_point(rng: &mut SmallRng) -> Point<EF> {
    Point::new((0..NUM_VARS).map(|_| rng.random::<EF>()).collect())
}

/// Independent multilinear evaluation `p(z) = Σ_x evals[x] · eq(x, z)`, big-endian
/// (coordinate 0 = most significant bit of x) — matching the PCS convention.
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
fn ed25519_tensor_pcs_roundtrip_and_matches_direct_eval() {
    let mut rng = SmallRng::seed_from_u64(1);
    let pcs = pcs();
    let evals = rand_evals(&mut rng);
    let point = rand_point(&mut rng);

    let mut ch_p = challenger();
    let (commitment, prover_data) = pcs.commit(evals.clone(), &mut ch_p);
    let (value, proof) = pcs.open(&prover_data, &point, &mut ch_p);

    assert_eq!(value, multilinear_eval(&evals, &point), "opened value wrong");

    let mut ch_v = challenger();
    pcs.verify(&commitment, &point, value, &proof, &mut ch_v)
        .expect("honest opening must verify over ed25519");
}

#[test]
fn ed25519_tensor_pcs_rejects_tampered_column() {
    let mut rng = SmallRng::seed_from_u64(3);
    let pcs = pcs();
    let evals = rand_evals(&mut rng);
    let point = rand_point(&mut rng);

    let mut ch_p = challenger();
    let (commitment, prover_data) = pcs.commit(evals, &mut ch_p);
    let (value, mut proof) = pcs.open(&prover_data, &point, &mut ch_p);

    // Corrupt one opened column value; its Merkle path no longer matches.
    proof.openings[0].opened_values[0][0] += F::ONE;

    let mut ch_v = challenger();
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut ch_v)
        .expect_err("a tampered column opening must be rejected");
    assert!(matches!(err, TensorPcsError::Merkle(_)), "got {err:?}");
}

#[test]
fn ed25519_tensor_pcs_rejects_wrong_value() {
    let mut rng = SmallRng::seed_from_u64(2);
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
    assert!(matches!(err, TensorPcsError::ClaimedValueMismatch), "got {err:?}");
}
