//! Quasar commit / open / verify over ed25519, at security-derived query
//! counts.
//!
//! Complements `pcs-profile` (single-shot, with exact op counters) with proper
//! criterion statistics. Run:
//! `cargo bench -p p3-quasar --bench quasar_pcs --features parallel`

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use p3_basefold::{BasefoldPcs, RandomFoldableCode};
use p3_challenger::DuplexChallengerPF;
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519, Poseidon2Ed25519,
};
use p3_merkle_tree::MerkleTreeMmcs;
use p3_qa_code::QaCode;
use p3_quasar::{QuasarParams, QuasarPcs};
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

/// `log₂N₁` — the paper's own benchmark shape.
const LOG_ROWS: usize = 6;
const BF_INV_RATE: usize = 4;
const BF_QUERIES: usize = 288;

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

fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("quasar_ed25519");
    group.sample_size(10);

    for num_vars in [14, 16, 18] {
        for index in [2, 4] {
            let log_cols = num_vars - LOG_ROWS;
            let num_queries = QuasarParams::from_security(LOG_ROWS, log_cols, index, 255, 103, 103)
                .expect("ed25519 admits a QA distance")
                .num_queries
                .min(index << log_cols);
            let params = QuasarParams::new(LOG_ROWS, log_cols, index, num_queries);
            let code = QaCode::<F>::new_random(
                log_cols,
                index,
                &mut ChaCha8Rng::seed_from_u64(0x0A_C0DE),
            );
            let basefold = BasefoldPcs::new(
                RandomFoldableCode::<F>::new(log_cols, BF_INV_RATE, 0xBA5E_F01D),
                mmcs(),
                BF_QUERIES.min(BF_INV_RATE << log_cols),
            );
            let (pcs, pk) = QuasarPcs::new(params, code, mmcs(), basefold);

            let mut rng = ChaCha8Rng::seed_from_u64(1);
            let evals: Vec<F> = (0..(1usize << num_vars)).map(|_| rng.random()).collect();
            let point: Vec<F> = (0..num_vars).map(|_| rng.random()).collect();

            let id = format!("rate_1_over_{index}");
            group.bench_with_input(BenchmarkId::new(format!("commit/{id}"), num_vars), &num_vars, |b, _| {
                b.iter(|| pcs.commit(evals.clone(), &mut challenger()))
            });

            let mut ch = challenger();
            let (commitment, pd) = pcs.commit(evals.clone(), &mut ch);
            let saved = ch.clone();
            group.bench_with_input(BenchmarkId::new(format!("open/{id}"), num_vars), &num_vars, |b, _| {
                b.iter(|| pcs.open(&pk, &pd, &point, &mut saved.clone()))
            });

            let (value, proof) = pcs.open(&pk, &pd, &point, &mut ch);
            group.bench_with_input(BenchmarkId::new(format!("verify/{id}"), num_vars), &num_vars, |b, _| {
                b.iter(|| {
                    pcs.verify(&commitment, &point, value, &proof, &mut challenger())
                        .expect("verify")
                })
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
