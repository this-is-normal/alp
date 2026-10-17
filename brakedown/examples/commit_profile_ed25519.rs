//! Single-run **stage profiler** for the Brakedown `commit` over ed25519.
//!
//! Unlike the criterion timing bench (which measures `commit`, `encode`, and
//! `merkle` as *separate* benchmarks — re-running work, with stage numbers from
//! different runs that only approximately sum), this profiles the stages **within
//! one commit execution**, so they sum exactly to that run's total and nothing is
//! executed twice. It's the PCS analogue of the multi-stark `tracing_forest` tree.
//!
//! `commit(evals)` = reshape/transpose + `code.encode` (linear-code encoding) +
//! `mmcs.commit_matrix` (Merkle/Poseidon2 hashing). We time each in place.
//!
//! Averaged over a few iterations for stability (per-stage averages still sum to
//! the average total). Run into a LOCAL target dir; add `--features parallel` to
//! see the breakdown with the Merkle build multithreaded:
//! ```text
//! CARGO_TARGET_DIR=/home/ubuntu/alp-target cargo run -p p3-brakedown \
//!   --example commit_profile_ed25519 --release [--features parallel]
//! ```

use std::time::{Duration, Instant};

use p3_brakedown::{BrakedownCode, LinearCode};
use p3_commit::Mmcs;
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519, Poseidon2Ed25519,
};
use p3_matrix::dense::RowMajorMatrix;
use p3_merkle_tree::MerkleTreeMmcs;
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

const CODE_SEED: u64 = 0xBADC0DE;
const PERM_SEED: u64 = 0xED25519_C0FFEE;

/// `(log_rows, log_cols)` shapes; square-ish. 2^12, 2^14, 2^16 cells.
const SHAPES: &[(usize, usize)] = &[(6, 6), (7, 7), (8, 8)];
const WARMUP: usize = 2;
const ITERS: usize = 10;

fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(PERM_SEED);
    Poseidon2Ed25519::<WIDTH>::new_from_rng(
        ED25519_POSEIDON2_FULL_ROUNDS,
        ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
        &mut rng,
    )
}

fn ms(d: Duration) -> f64 {
    d.as_secs_f64() * 1e3
}

fn main() {
    println!("Brakedown commit stage profiler over ed25519 (single-run decomposition; avg of {ITERS} iters)");
    #[cfg(feature = "parallel")]
    println!("parallelism: ON (rayon via maybe-rayon)");
    #[cfg(not(feature = "parallel"))]
    println!("parallelism: OFF (build with --features parallel to enable)");
    println!(
        "\n{:>4} {:>6} {:>6} {:>7} | {:>10} {:>10} {:>10} | {:>10}",
        "vars", "rows", "cols", "code_n", "reshape", "encode", "merkle", "commit"
    );

    for &(log_rows, log_cols) in SHAPES {
        let (rows, cols) = (1usize << log_rows, 1usize << log_cols);
        let num_vars = log_rows + log_cols;

        let mut code_rng = SmallRng::seed_from_u64(CODE_SEED);
        let code = BrakedownCode::<F, EF>::from_message_len(cols, &mut code_rng);
        let code_n = code.codeword_len();
        let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);

        let mut rng = SmallRng::seed_from_u64(1);
        let make_evals = |rng: &mut SmallRng| -> Vec<F> {
            (0..(rows * cols)).map(|_| rng.random::<F>()).collect()
        };

        let (mut t_reshape, mut t_encode, mut t_merkle) =
            (Duration::ZERO, Duration::ZERO, Duration::ZERO);

        for i in 0..(WARMUP + ITERS) {
            let evals = make_evals(&mut rng);

            // ---- the three stages of one commit, timed in place ----
            let t0 = Instant::now();
            let m_natural = RowMajorMatrix::new(evals, cols);
            let msg = m_natural.transpose(); // cols × rows
            let d_reshape = t0.elapsed();

            let t1 = Instant::now();
            let encoded = code.encode(&msg); // n × rows
            let d_encode = t1.elapsed();

            let t2 = Instant::now();
            let _committed = mmcs.commit_matrix(encoded); // Merkle build (hashing)
            let d_merkle = t2.elapsed();

            if i >= WARMUP {
                t_reshape += d_reshape;
                t_encode += d_encode;
                t_merkle += d_merkle;
            }
        }

        let n = ITERS as u32;
        let (r, e, m) = (t_reshape / n, t_encode / n, t_merkle / n);
        let total = r + e + m;
        println!(
            "{:>4} {:>6} {:>6} {:>7} | {:>9.3}m {:>9.3}m {:>9.3}m | {:>9.3}m",
            num_vars,
            rows,
            cols,
            code_n,
            ms(r),
            ms(e),
            ms(m),
            ms(total),
        );
    }

    println!("\n(reshape + encode + merkle sum exactly to commit — one execution, no re-runs)");
    println!("merkle = Poseidon2 Merkle build (the ~97% hashing); encode = Brakedown linear code");
}
