//! Criterion timing for the Brakedown tensor PCS **over the native ed25519 base
//! field** (`p = 2^255 − 19`): `commit`, `open`, `verify` measured separately
//! across polynomial sizes. ed25519 counterpart to `tensor_pcs.rs` (BabyBear).
//!
//! Uses the plain (un-instrumented) `TensorPcs`, so no counting overhead pollutes
//! the wall-clock. Pair with the exact op counts from
//! `examples/pcs_cost_sweep_ed25519.rs` to drive the cost model on the real field.
//!
//! Run (into a LOCAL target dir — sshfs can't link):
//! ```text
//! CARGO_TARGET_DIR=/home/ubuntu/alp-target \
//!   cargo bench -p p3-brakedown --bench tensor_pcs_ed25519
//! ```

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use p3_brakedown::{BrakedownCode, TensorPcs};
use p3_challenger::DuplexChallengerPF;
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519, Poseidon2Ed25519,
};
use p3_merkle_tree::MerkleTreeMmcs;
use p3_multilinear_util::point::Point;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

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
type Pcs = TensorPcs<F, EF, BrakedownCode<F, EF>, MyMmcs>;

const CODE_SEED: u64 = 0xBADC0DE;
const PERM_SEED: u64 = 0xED25519_C0FFEE;
const NUM_QUERIES: usize = 32;

/// `(log_rows, log_cols)` shapes; square-ish so encoding and query cost balance.
/// ed25519 muls are ~10-35x BabyBear's, so we stop at 2^16 to keep the bench brisk.
const SHAPES: &[(usize, usize)] = &[(6, 6), (7, 7), (8, 8)]; // 2^12, 2^14, 2^16 cells

fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(PERM_SEED);
    Poseidon2Ed25519::<WIDTH>::new_from_rng(
        ED25519_POSEIDON2_FULL_ROUNDS,
        ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
        &mut rng,
    )
}

fn build_pcs(log_rows: usize, cols: usize) -> Pcs {
    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);
    let mut code_rng = SmallRng::seed_from_u64(CODE_SEED);
    let code = BrakedownCode::<F, EF>::from_message_len(cols, &mut code_rng);
    TensorPcs::new(code, mmcs, log_rows, NUM_QUERIES)
}

fn rand_evals(num_vars: usize) -> Vec<F> {
    let mut rng = SmallRng::seed_from_u64(1);
    (0..(1usize << num_vars)).map(|_| rng.random::<F>()).collect()
}

fn rand_point(num_vars: usize) -> Point<EF> {
    let mut rng = SmallRng::seed_from_u64(2);
    Point::new((0..num_vars).map(|_| rng.random::<EF>()).collect())
}

fn bench_commit(c: &mut Criterion) {
    let mut g = c.benchmark_group("tensor_pcs_ed25519/commit");
    for &(log_rows, log_cols) in SHAPES {
        let (cols, num_vars) = (1usize << log_cols, log_rows + log_cols);
        let pcs = build_pcs(log_rows, cols);
        let evals = rand_evals(num_vars);
        g.bench_with_input(BenchmarkId::from_parameter(num_vars), &num_vars, |b, _| {
            b.iter_batched(
                || evals.clone(),
                |evals| pcs.commit(evals, &mut Challenger::new(perm())),
                BatchSize::SmallInput,
            );
        });
    }
    g.finish();
}

fn bench_open(c: &mut Criterion) {
    let mut g = c.benchmark_group("tensor_pcs_ed25519/open");
    for &(log_rows, log_cols) in SHAPES {
        let (cols, num_vars) = (1usize << log_cols, log_rows + log_cols);
        let pcs = build_pcs(log_rows, cols);
        let point = rand_point(num_vars);
        let (_, prover_data) = pcs.commit(rand_evals(num_vars), &mut Challenger::new(perm()));
        g.bench_with_input(BenchmarkId::from_parameter(num_vars), &num_vars, |b, _| {
            b.iter(|| pcs.open(&prover_data, &point, &mut Challenger::new(perm())));
        });
    }
    g.finish();
}

fn bench_verify(c: &mut Criterion) {
    let mut g = c.benchmark_group("tensor_pcs_ed25519/verify");
    for &(log_rows, log_cols) in SHAPES {
        let (cols, num_vars) = (1usize << log_cols, log_rows + log_cols);
        let pcs = build_pcs(log_rows, cols);
        let point = rand_point(num_vars);
        let mut ch_p = Challenger::new(perm());
        let (commitment, prover_data) = pcs.commit(rand_evals(num_vars), &mut ch_p);
        let (value, proof) = pcs.open(&prover_data, &point, &mut ch_p);
        g.bench_with_input(BenchmarkId::from_parameter(num_vars), &num_vars, |b, _| {
            b.iter(|| {
                pcs.verify(&commitment, &point, value, &proof, &mut Challenger::new(perm()))
                    .expect("honest opening must verify");
            });
        });
    }
    g.finish();
}

// The `commit` stage split (encode vs Merkle hashing) is measured by the
// single-run `commit_profile_ed25519` example, which decomposes ONE commit so
// the stages sum exactly — rather than re-running each stage as its own bench.

criterion_group!(tensor_pcs_ed25519, bench_commit, bench_open, bench_verify);
criterion_main!(tensor_pcs_ed25519);
