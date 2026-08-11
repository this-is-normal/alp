//! Criterion timing for the Brakedown tensor PCS: `commit`, `open`, `verify`
//! measured **separately** across a few polynomial sizes.
//!
//! This is the *timed* counterpart to `examples/pcs_cost_sweep.rs`: it uses the
//! plain (un-instrumented) `TensorPcs` so no counting overhead pollutes the
//! wall-clock numbers. Pair the per-op times here with the deterministic op
//! counts from the sweep to drive the cost model.
//!
//! Run (into a LOCAL target dir — sshfs can't link):
//! ```text
//! CARGO_TARGET_DIR=/home/ubuntu/alp-target \
//!   cargo bench -p p3-brakedown --bench tensor_pcs
//! ```
//! Filter one stage:  `... --bench tensor_pcs -- commit`

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_brakedown::{BrakedownCode, TensorPcs};
use p3_challenger::DuplexChallenger;
use p3_field::extension::BinomialExtensionField;
use p3_field::Field;
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

const CODE_SEED: u64 = 0xBADC0DE;
const PERM_SEED: u64 = 0xC0FFEE;
const NUM_QUERIES: usize = 32;

/// `(log_rows, log_cols)` shapes; square-ish so encoding and query cost are balanced.
const SHAPES: &[(usize, usize)] = &[(6, 6), (8, 8), (9, 9)]; // 2^12, 2^16, 2^18 cells

fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(PERM_SEED);
    Perm::new_from_rng_128(&mut rng)
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
    let mut g = c.benchmark_group("tensor_pcs/commit");
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
    let mut g = c.benchmark_group("tensor_pcs/open");
    for &(log_rows, log_cols) in SHAPES {
        let (cols, num_vars) = (1usize << log_cols, log_rows + log_cols);
        let pcs = build_pcs(log_rows, cols);
        let point = rand_point(num_vars);
        // Commit once; each iteration re-opens the same committed data.
        let (_, prover_data) = pcs.commit(rand_evals(num_vars), &mut Challenger::new(perm()));
        g.bench_with_input(BenchmarkId::from_parameter(num_vars), &num_vars, |b, _| {
            b.iter(|| pcs.open(&prover_data, &point, &mut Challenger::new(perm())));
        });
    }
    g.finish();
}

fn bench_verify(c: &mut Criterion) {
    let mut g = c.benchmark_group("tensor_pcs/verify");
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

criterion_group!(tensor_pcs, bench_commit, bench_open, bench_verify);
criterion_main!(tensor_pcs);
