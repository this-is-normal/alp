//! End-to-end Quasar over the native ed25519 base field (`p = 2²⁵⁵ − 19`):
//! commit, open, verify, and rejection of every part of the proof a cheating
//! prover would have to forge.

use p3_basefold::{BasefoldPcs, RandomFoldableCode};
use p3_challenger::DuplexChallengerPF;
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519, Poseidon2Ed25519,
};
use p3_field::PrimeCharacteristicRing;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_qa_code::QaCode;
use p3_quasar::{QuasarError, QuasarParams, QuasarPcs, QuasarProverKey, eval_mle};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;

type F = Ed25519;

const WIDTH: usize = 3;
const RATE: usize = 2;
const OUT: usize = 1;

type Perm = Poseidon2Ed25519<WIDTH>;
type MyHash = PaddingFreeSponge<Perm, WIDTH, RATE, OUT>;
type MyCompress = TruncatedPermutation<Perm, 2, OUT, WIDTH>;
type MyMmcs = MerkleTreeMmcs<Ed25519, Ed25519, MyHash, MyCompress, 2, OUT>;
type Challenger = DuplexChallengerPF<Ed25519, Perm, WIDTH, RATE>;
type Pcs = QuasarPcs<F, MyMmcs>;

/// Column queries. A real instance opens ~425 columns (see
/// `column_count_for_a_real_instance`); these test codewords are far too short
/// for that, so the tests use a token count.
const NUM_QUERIES: usize = 16;
/// BaseFold query repetitions for the code-switching openings.
const BF_QUERIES: usize = 12;
/// BaseFold inverse rate.
const BF_INV_RATE: usize = 4;

fn perm() -> Perm {
    let mut rng = ChaCha8Rng::seed_from_u64(0xED25519_C0FFEE);
    Poseidon2Ed25519::<WIDTH>::new_from_rng(
        ED25519_POSEIDON2_FULL_ROUNDS,
        ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
        &mut rng,
    )
}

fn mmcs() -> MyMmcs {
    let p = perm();
    MyMmcs::new(MyHash::new(p.clone()), MyCompress::new(p), 0)
}

fn challenger() -> Challenger {
    Challenger::new(perm())
}

fn build(log_rows: usize, log_cols: usize, index: usize) -> (Pcs, QuasarProverKey<F, MyMmcs>) {
    let params = QuasarParams::new(log_rows, log_cols, index, NUM_QUERIES);
    let code = QaCode::<F>::new_random(log_cols, index, &mut ChaCha8Rng::seed_from_u64(0x0A_C0DE));
    let basefold = BasefoldPcs::new(
        RandomFoldableCode::<F>::new(log_cols, BF_INV_RATE, 0xBA5E_F01D),
        mmcs(),
        BF_QUERIES,
    );
    QuasarPcs::new(params, code, mmcs(), basefold)
}

fn rand_evals(rng: &mut ChaCha8Rng, num_vars: usize) -> Vec<F> {
    (0..(1 << num_vars)).map(|_| rng.random::<F>()).collect()
}

fn rand_point(rng: &mut ChaCha8Rng, num_vars: usize) -> Vec<F> {
    (0..num_vars).map(|_| rng.random::<F>()).collect()
}

/// Commit / open / verify at several shapes and both rates, checking the opened
/// value against an independent multilinear evaluation.
#[test]
fn roundtrip_and_matches_direct_eval() {
    for (log_rows, log_cols, index) in [(3, 5, 2), (2, 6, 2), (4, 4, 2), (3, 5, 4), (1, 6, 4)] {
        let (pcs, pk) = build(log_rows, log_cols, index);
        let num_vars = log_rows + log_cols;
        let mut rng = ChaCha8Rng::seed_from_u64(1 + num_vars as u64 + index as u64);
        let evals = rand_evals(&mut rng, num_vars);
        let point = rand_point(&mut rng, num_vars);

        let mut ch = challenger();
        let (commitment, pd) = pcs.commit(evals.clone(), &mut ch);
        let (value, proof) = pcs.open(&pk, &pd, &point, &mut ch);

        assert_eq!(
            value,
            eval_mle(&evals, &point),
            "opened value wrong at {log_rows}x{log_cols} index {index}"
        );

        pcs.verify(&commitment, &point, value, &proof, &mut challenger())
            .unwrap_or_else(|e| {
                panic!("honest reject at {log_rows}x{log_cols} index {index}: {e:?}")
            });
    }
}

/// The one-row degenerate shape: `N₁ = 1` makes the proximity row equal a
/// scalar multiple of the single matrix row.
#[test]
fn roundtrip_with_a_single_row() {
    let (pcs, pk) = build(0, 6, 2);
    let mut rng = ChaCha8Rng::seed_from_u64(11);
    let evals = rand_evals(&mut rng, 6);
    let point = rand_point(&mut rng, 6);

    let mut ch = challenger();
    let (commitment, pd) = pcs.commit(evals.clone(), &mut ch);
    let (value, proof) = pcs.open(&pk, &pd, &point, &mut ch);
    assert_eq!(value, eval_mle(&evals, &point));
    pcs.verify(&commitment, &point, value, &proof, &mut challenger())
        .expect("single-row instance must verify");
}

fn honest_proof(
    seed: u64,
) -> (
    Pcs,
    <MyMmcs as p3_commit::Mmcs<F>>::Commitment,
    Vec<F>,
    F,
    p3_quasar::QuasarProof<F, MyMmcs>,
) {
    let (pcs, pk) = build(3, 5, 2);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let evals = rand_evals(&mut rng, 8);
    let point = rand_point(&mut rng, 8);
    let mut ch = challenger();
    let (commitment, pd) = pcs.commit(evals, &mut ch);
    let (value, proof) = pcs.open(&pk, &pd, &point, &mut ch);
    (pcs, commitment, point, value, proof)
}

#[test]
fn rejects_wrong_claimed_value() {
    let (pcs, commitment, point, value, proof) = honest_proof(2);
    let err = pcs
        .verify(
            &commitment,
            &point,
            value + F::ONE,
            &proof,
            &mut challenger(),
        )
        .expect_err("a wrong claimed value must be rejected");
    assert!(matches!(err, QuasarError::ClaimedValueMismatch), "got {err:?}");
}

/// Forging the value *inside* the proof does not help: the claimed value only
/// enters through the final `q₁(z_R)` claim, so the reduction sumcheck catches
/// it.
#[test]
fn rejects_forged_claimed_value_in_proof() {
    let (pcs, commitment, point, value, mut proof) = honest_proof(12);
    proof.claimed_value += F::ONE;
    let err = pcs
        .verify(&commitment, &point, value + F::ONE, &proof, &mut challenger())
        .expect_err("a forged claimed value must be rejected");
    assert!(matches!(err, QuasarError::ReductionClaim), "got {err:?}");
}

#[test]
fn rejects_tampered_column() {
    let (pcs, commitment, point, value, mut proof) = honest_proof(3);
    proof.column_openings[0][0][0] += F::ONE;
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut challenger())
        .expect_err("a tampered column must be rejected");
    assert!(matches!(err, QuasarError::Merkle(_)), "got {err:?}");
}

/// Corrupting a block evaluation breaks the selector sumcheck's terminal check:
/// this is the link between the Merkle-committed columns and the BaseFold-
/// committed blocks.
#[test]
fn rejects_tampered_block_evaluation() {
    let (pcs, commitment, point, value, mut proof) = honest_proof(4);
    proof.block_evals_at_e1[0] += F::ONE;
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut challenger())
        .expect_err("a tampered block evaluation must be rejected");
    assert!(matches!(err, QuasarError::SelectorTerminal), "got {err:?}");
}

/// Corrupting a parity-block evaluation breaks the QA-encoding claim — the
/// code-switching check proper.
#[test]
fn rejects_tampered_parity_evaluation() {
    let (pcs, commitment, point, value, mut proof) = honest_proof(5);
    proof.parity_evals_at_gamma[0] += F::ONE;
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut challenger())
        .expect_err("a tampered parity evaluation must be rejected");
    assert!(matches!(err, QuasarError::EncodingClaim), "got {err:?}");
}

/// Corrupting the preprocessed coefficient evaluation breaks the outer
/// encoding sumcheck's terminal check.
#[test]
fn rejects_tampered_coefficient_evaluation() {
    let (pcs, commitment, point, value, mut proof) = honest_proof(6);
    proof.coeff_evals_at_e2[0] += F::ONE;
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut challenger())
        .expect_err("a tampered coefficient evaluation must be rejected");
    assert!(
        matches!(err, QuasarError::EncodingOuterTerminal),
        "got {err:?}"
    );
}

/// Corrupting the encoding witness's reduced evaluation breaks the inner layer.
#[test]
fn rejects_tampered_input_evaluation() {
    let (pcs, commitment, point, value, mut proof) = honest_proof(7);
    proof.input_evals_at_e3[0] += F::ONE;
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut challenger())
        .expect_err("a tampered input evaluation must be rejected");
    assert!(
        matches!(err, QuasarError::EncodingInnerTerminal),
        "got {err:?}"
    );
}

/// Corrupting a value at the reduction point breaks the batched BaseFold
/// opening: the final tie between all the claims and the committed blocks.
#[test]
fn rejects_tampered_evaluation_at_zeta() {
    let (pcs, commitment, point, value, mut proof) = honest_proof(8);
    proof.evals_at_zeta[0] += F::ONE;
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut challenger())
        .expect_err("a tampered zeta evaluation must be rejected");
    assert!(
        matches!(
            err,
            QuasarError::ReductionTerminal | QuasarError::OpeningValue
        ),
        "got {err:?}"
    );
}

/// Corrupting a BaseFold leaf must be caught by the batched opening.
#[test]
fn rejects_tampered_basefold_leaf() {
    let (pcs, commitment, point, value, mut proof) = honest_proof(9);
    proof.online_opening.queries[0].opened_values[0][0][0] += F::ONE;
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut challenger())
        .expect_err("a tampered BaseFold leaf must be rejected");
    assert!(matches!(err, QuasarError::Basefold(_)), "got {err:?}");
}

/// A proof for one point must not verify at another.
#[test]
fn rejects_wrong_point() {
    let (pcs, commitment, mut point, value, proof) = honest_proof(10);
    point[0] += F::ONE;
    assert!(
        pcs.verify(&commitment, &point, value, &proof, &mut challenger())
            .is_err(),
        "a proof must not transfer to a different point"
    );
}

/// What a production instance over ed25519 actually opens, from the proven QA
/// distance — recorded here because it is what makes proof size depend so
/// strongly on `N₁`.
#[test]
fn column_count_for_a_real_instance() {
    let p = QuasarParams::from_security(6, 16, 2, 255, 103, 103).unwrap();
    assert_eq!(p.num_queries, 425);
    let p4 = QuasarParams::from_security(6, 16, 4, 255, 103, 103).unwrap();
    assert_eq!(p4.num_queries, 256);
}

/// Grinding must gate the column sample: a ground proof verifies, a proof with
/// the wrong witness does not.
#[test]
fn grinding_gates_the_column_sample() {
    const POW: usize = 8;
    let params = QuasarParams::new(3, 5, 2, NUM_QUERIES);
    let code = QaCode::<F>::new_random(5, 2, &mut ChaCha8Rng::seed_from_u64(0x0A_C0DE));
    let basefold = BasefoldPcs::new_with_grinding(
        RandomFoldableCode::<F>::new(5, BF_INV_RATE, 0xBA5E_F01D),
        mmcs(),
        BF_QUERIES.min(BF_INV_RATE << 5),
        POW,
    );
    let (pcs, pk) = QuasarPcs::new(params, code, mmcs(), basefold);
    let pcs = pcs.with_grinding(POW);

    let mut rng = ChaCha8Rng::seed_from_u64(30);
    let evals = rand_evals(&mut rng, 8);
    let point = rand_point(&mut rng, 8);

    let mut ch = challenger();
    let (commitment, pd) = pcs.commit(evals.clone(), &mut ch);
    let (value, mut proof) = pcs.open(&pk, &pd, &point, &mut ch);
    assert_eq!(value, eval_mle(&evals, &point));
    pcs.verify(&commitment, &point, value, &proof, &mut challenger())
        .expect("a ground proof must verify");

    proof.pow_witness += F::ONE;
    let err = pcs
        .verify(&commitment, &point, value, &proof, &mut challenger())
        .expect_err("a bad proof-of-work witness must be rejected");
    assert!(matches!(err, QuasarError::InvalidProofOfWork), "got {err:?}");
}

/// A verifier-only instance — built from the coefficient *commitment* rather
/// than the QA code, so it never runs the setup encode-and-commit — accepts
/// exactly what the full instance accepts, and rejects the same tampering.
///
/// This is what a recursion circuit builds: reconstructing the code there costs
/// `O(N₂)` field work and hashing per proof, and the verifier reads only the digest.
#[test]
fn verifier_only_instance_matches_the_full_one() {
    let (log_rows, log_cols, index) = (3, 5, 2);
    let (pcs, pk) = build(log_rows, log_cols, index);
    let num_vars = log_rows + log_cols;
    let mut rng = ChaCha8Rng::seed_from_u64(4242);
    let evals = rand_evals(&mut rng, num_vars);
    let point = rand_point(&mut rng, num_vars);

    let mut ch = challenger();
    let (commitment, pd) = pcs.commit(evals, &mut ch);
    let (value, proof) = pcs.open(&pk, &pd, &point, &mut ch);

    let verifier = QuasarPcs::new_verifier(
        QuasarParams::new(log_rows, log_cols, index, NUM_QUERIES),
        mmcs(),
        BasefoldPcs::new(
            RandomFoldableCode::<F>::new(log_cols, BF_INV_RATE, 0xBA5E_F01D),
            mmcs(),
            BF_QUERIES,
        ),
        pcs.coeff_commitment().clone(),
    );

    verifier
        .verify(&commitment, &point, value, &proof, &mut challenger())
        .expect("verifier-only instance must accept an honest proof");

    // ... and it is not accepting everything.
    verifier
        .verify(&commitment, &point, value + F::ONE, &proof, &mut challenger())
        .expect_err("verifier-only instance must reject a wrong value");
}

/// Why is the *coefficient* Basefold opening 3.4x the *online* one in the
/// OpenVM guest (42.9 M vs 12.5 M instructions at 16 signatures), when it
/// batches fewer polynomials at the same point with the same query count?
///
/// Prints the shape of both openings so the asymmetry can be attributed rather
/// than guessed at. Run with `--nocapture`.
#[test]
fn opening_shapes_online_vs_coefficient() {
    let (log_rows, log_cols, index) = (3, 6, 2);
    let (pcs, pk) = build(log_rows, log_cols, index);
    let num_vars = log_rows + log_cols;
    let mut rng = ChaCha8Rng::seed_from_u64(77);
    let evals = rand_evals(&mut rng, num_vars);
    let point = rand_point(&mut rng, num_vars);

    let mut ch = challenger();
    let (_commitment, pd) = pcs.commit(evals, &mut ch);
    let (_value, proof) = pcs.open(&pk, &pd, &point, &mut ch);

    let describe = |name: &str, opening: &p3_basefold::BasefoldProof<F, MyMmcs>| {
        let rounds = opening.queries.len();
        let mut total_opened = 0usize;
        let mut widths = Vec::new();
        for round in &opening.queries {
            let rows = round.opened_values.len();
            let width = round.opened_values.first().map_or(0, |v| v[0].len());
            total_opened += rows * width;
            widths.push((rows, width));
        }
        println!(
            "{name}: rounds={rounds} fold_roots={} sumcheck={} opened_elements={total_opened}",
            opening.fold_roots.len(),
            opening.sumcheck.len(),
        );
        println!("  (rows, width) per round: {widths:?}");
    };

    describe("online   ", &proof.online_opening);
    describe("coeff    ", &proof.coeff_opening);
}
