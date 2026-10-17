//! Tensor-PCS tests using the **real recursive Brakedown code** (M4).
//!
//! Mirrors `tensor_pcs.rs`, but instantiates `TensorPcs` with a genuine
//! `BrakedownCode::from_message_len` linear-time encoder instead of Reed–Solomon.
//! This makes the commitment truly Brakedown (linear-time encoding), not Ligero.
//!
//! Validates the real cryptographic behaviour:
//!  * honest roundtrip verifies AND the opened value equals an independently
//!    computed multilinear evaluation;
//!  * a tampered column opening is rejected.
//!
//! Note: for Brakedown the codeword length `n` is *not* a power of two — the
//! tensor PCS and Merkle tree handle arbitrary `n`; we only need `num_queries <= n`.

use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_brakedown::{BrakedownCode, LinearCode, TensorPcs, TensorPcsError};
use p3_challenger::DuplexChallenger;
use p3_field::extension::BinomialExtensionField;
use p3_field::{Field, PrimeCharacteristicRing};
use p3_merkle_tree::MerkleTreeMmcs;
use p3_multilinear_util::point::Point;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

type F = BabyBear;
type EF = BinomialExtensionField<F, 4>;
type Perm = Poseidon2BabyBear<16>;
type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;
type PackedF = <F as Field>::Packing;
type MyMmcs = MerkleTreeMmcs<PackedF, PackedF, MyHash, MyCompress, 2, 8>;
type Challenger = DuplexChallenger<F, Perm, 16, 8>;
type Pcs = TensorPcs<F, EF, BrakedownCode<F, EF>, MyMmcs>;

const LOG_ROWS: usize = 3; // rows = 8
const COLS: usize = 32; // Brakedown message length (power of two; recurses one level)
const NUM_QUERIES: usize = 12;
const NUM_VARS: usize = 8; // rows(8) * cols(32) = 256 = 2^8

fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(0xC0FFEE);
    Perm::new_from_rng_128(&mut rng)
}

fn challenger() -> Challenger {
    Challenger::new(perm())
}

/// Build the PCS with a real Brakedown code. The code's random A/B matrices are
/// drawn from a *fixed-seed* RNG so the (single) `Pcs` value used by both prover
/// and verifier shares one code.
fn pcs() -> Pcs {
    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);
    let mut code_rng = SmallRng::seed_from_u64(0xBADC0DE);
    let code = BrakedownCode::<F, EF>::from_message_len(COLS, &mut code_rng);
    assert_eq!(code.message_len(), COLS);
    // Sanity: brakedown codeword length is not a power of two, and we must be
    // able to sample the requested queries within it.
    assert!(NUM_QUERIES <= code.codeword_len());
    TensorPcs::new(code, mmcs, LOG_ROWS, NUM_QUERIES)
}

fn rand_evals(rng: &mut SmallRng) -> Vec<F> {
    (0..(1 << NUM_VARS)).map(|_| rng.random::<F>()).collect()
}

fn rand_point(rng: &mut SmallRng) -> Point<EF> {
    Point::new((0..NUM_VARS).map(|_| rng.random::<EF>()).collect())
}

/// Independent multilinear evaluation `p(z) = Σ_x evals[x] · eq(x, z)`,
/// big-endian (coordinate 0 = most significant bit of x) — matching the PCS.
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
fn brakedown_tensor_pcs_roundtrip_and_matches_direct_eval() {
    let mut rng = SmallRng::seed_from_u64(1);
    let pcs = pcs();
    let evals = rand_evals(&mut rng);
    let point = rand_point(&mut rng);

    // Prover: commit then open (one challenger).
    let mut ch_p = challenger();
    let (commitment, prover_data) = pcs.commit(evals.clone(), &mut ch_p);
    let (value, proof) = pcs.open(&prover_data, &point, &mut ch_p);

    // The opened value must equal the true multilinear evaluation.
    assert_eq!(value, multilinear_eval(&evals, &point), "opened value wrong");

    // Verifier: fresh identically-seeded challenger.
    let mut ch_v = challenger();
    pcs.verify(&commitment, &point, value, &proof, &mut ch_v)
        .expect("honest opening must verify under the real Brakedown code");
}

#[test]
fn brakedown_tensor_pcs_rejects_tampered_column() {
    let mut rng = SmallRng::seed_from_u64(3);
    let pcs = pcs();
    let evals = rand_evals(&mut rng);
    let point = rand_point(&mut rng);

    let mut ch_p = challenger();
    let (commitment, prover_data) = pcs.commit(evals, &mut ch_p);
    let (value, mut proof) = pcs.open(&prover_data, &point, &mut ch_p);

    // Mutation: corrupt one opened column value; its Merkle path no longer matches.
    proof.openings[0][0][0] += F::ONE;

    let mut ch_v = challenger();
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut ch_v)
        .expect_err("a tampered column opening must be rejected");
    // Merkle verification runs before the algebraic checks.
    assert!(matches!(err, TensorPcsError::Merkle(_)), "got {err:?}");
}

#[test]
fn brakedown_tensor_pcs_rejects_wrong_value() {
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
