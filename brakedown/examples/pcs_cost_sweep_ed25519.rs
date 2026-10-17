//! Brakedown/Ligero tensor-PCS **cost-model sweep over the native ed25519 base
//! field** (`p = 2^255 − 19`). ed25519 counterpart to `pcs_cost_sweep.rs`.
//!
//! Drives the tensor `TensorPcs` through an *instrumented* hash stack and prints,
//! per configuration, the prover/verifier **permutation counts** plus serialized
//! proof size — the deterministic `op_count` half of `time ≈ Σ op_count · per_op_time`.
//!
//! We count the **underlying permutation** (`CountingPerm`), not hasher calls, so
//! a wide-column leaf hash is charged its true `⌈rows/RATE⌉` permutations. Roles:
//!   * **L** = leaf-sponge permutations, **C** = compression permutations,
//!   * **S** = challenger-sponge permutations (Fiat-Shamir).
//! Their sum is the total permutation count. ed25519 has `Packing = Self`, so
//! each count is an **exact scalar permutation count**.
//!
//! The final section is a **Brakedown vs Ligero (Reed–Solomon)** head-to-head at
//! matching shape/queries — same measurement harness, only the linear code differs.
//!
//! Run (into a LOCAL target dir — sshfs can't link):
//! ```text
//! CARGO_TARGET_DIR=/home/ubuntu/alp-target \
//!   cargo run -p p3-brakedown --example pcs_cost_sweep_ed25519 --release
//! ```

use p3_brakedown::{BrakedownCode, LinearCodeExt, ReedSolomonCode, TensorPcs};
use p3_challenger::DuplexChallengerPF;
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519, Poseidon2Ed25519,
};
use p3_merkle_tree::MerkleTreeMmcs;
use p3_multilinear_util::point::Point;
use p3_stark_stats::{serialized_size, CountingPerm, HashCounts};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

type F = Ed25519;
type EF = Ed25519;

const WIDTH: usize = 3;
const RATE: usize = 2;
const OUT: usize = 1;

type Perm = Poseidon2Ed25519<WIDTH>;
type CPerm = CountingPerm<Perm>;
type MyHash = PaddingFreeSponge<CPerm, WIDTH, RATE, OUT>;
type MyCompress = TruncatedPermutation<CPerm, 2, OUT, WIDTH>;
type MyMmcs = MerkleTreeMmcs<Ed25519, Ed25519, MyHash, MyCompress, 2, OUT>;
type MyChallenger = DuplexChallengerPF<Ed25519, CPerm, WIDTH, RATE>;

const CODE_SEED: u64 = 0xBADC0DE;
const PERM_SEED: u64 = 0xED25519_C0FFEE;
/// Reed–Solomon inverse rate `n/k` (distance δ = 1 − 1/rate).
const RS_RATE: usize = 2;

fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(PERM_SEED);
    Poseidon2Ed25519::<WIDTH>::new_from_rng(
        ED25519_POSEIDON2_FULL_ROUNDS,
        ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
        &mut rng,
    )
}

fn counting_challenger(counts: &HashCounts) -> MyChallenger {
    DuplexChallengerPF::new(CountingPerm::with_counter(perm(), counts.sponge.clone()))
}

struct Row {
    label: &'static str,
    num_vars: usize,
    rows: usize,
    cols: usize,
    code_n: usize,
    queries: usize,
    prover: (usize, usize, usize),
    verifier: (usize, usize, usize),
    proof_bytes: usize,
}

fn bd_code(cols: usize) -> BrakedownCode<F, EF> {
    let mut r = SmallRng::seed_from_u64(CODE_SEED);
    BrakedownCode::<F, EF>::from_message_len(cols, &mut r)
}

fn rs_code(cols: usize) -> ReedSolomonCode<F> {
    ReedSolomonCode::new(cols, RS_RATE * cols)
}

/// Measure one PCS configuration — generic over the linear code, so Brakedown and
/// Reed–Solomon go through the identical counting harness. `requested_queries` is
/// clamped to the codeword length.
fn measure_code<C: LinearCodeExt<F, EF>>(
    label: &'static str,
    code: C,
    log_rows: usize,
    requested_queries: usize,
) -> Row {
    let cols = code.message_len();
    let code_n = code.codeword_len();
    let num_vars = log_rows + cols.trailing_zeros() as usize;
    let num_queries = requested_queries.min(code_n);

    let counts = HashCounts::new();
    let hash = MyHash::new(CountingPerm::with_counter(perm(), counts.leaf_hash.clone()));
    let compress = MyCompress::new(CountingPerm::with_counter(perm(), counts.compress.clone()));
    let mmcs = MyMmcs::new(hash, compress, 0);
    let pcs = TensorPcs::<F, EF, C, MyMmcs>::new(code, mmcs, log_rows, num_queries);

    let mut rng = SmallRng::seed_from_u64(1);
    let evals: Vec<F> = (0..(1usize << num_vars)).map(|_| rng.random::<F>()).collect();
    let point = Point::new((0..num_vars).map(|_| rng.random::<EF>()).collect());

    counts.reset();
    let mut ch_p = counting_challenger(&counts);
    let (commitment, prover_data) = pcs.commit(evals, &mut ch_p);
    let (value, proof) = pcs.open(&prover_data, &point, &mut ch_p);
    let prover = counts.snapshot();
    let proof_bytes = serialized_size(&proof);

    counts.reset();
    let mut ch_v = counting_challenger(&counts);
    pcs.verify(&commitment, &point, value, &proof, &mut ch_v)
        .expect("honest opening must verify");
    let verifier = counts.snapshot();

    Row {
        label,
        num_vars,
        rows: 1 << log_rows,
        cols,
        code_n,
        queries: num_queries,
        prover,
        verifier,
        proof_bytes,
    }
}

/// Brakedown at a square split with `q=32`.
fn measure(log_rows: usize, cols: usize, num_queries: usize) -> Row {
    measure_code("brakedown", bd_code(cols), log_rows, num_queries)
}

fn header(title: &str) {
    println!("\n=== {title} ===");
    println!(
        "{:>10} {:>4} {:>6} {:>6} {:>7} {:>4} | {:>22} | {:>22} | {:>9}",
        "code", "vars", "rows", "cols", "code_n", "q", "prover perms L/C/S", "verifier perms L/C/S",
        "proof(B)"
    );
}

fn fmt_perms((l, c, s): (usize, usize, usize)) -> String {
    format!("{l}/{c}/{s} ={}", l + c + s)
}

fn print_row(r: &Row) {
    println!(
        "{:>10} {:>4} {:>6} {:>6} {:>7} {:>4} | {:>22} | {:>22} | {:>9}",
        r.label,
        r.num_vars,
        r.rows,
        r.cols,
        r.code_n,
        r.queries,
        fmt_perms(r.prover),
        fmt_perms(r.verifier),
        r.proof_bytes,
    );
}

fn main() {
    println!("Brakedown/Ligero tensor-PCS permutation-cost sweep over ed25519 (Poseidon2 + Merkle)");
    println!("counts are POSEIDON2 PERMUTATIONS by role: L=leaf-sponge  C=compress  S=challenger (=sum total)");

    header("scale: num_vars up, square split, q=32 (Brakedown)");
    for &lv in &[8usize, 12, 14, 16, 18] {
        let log_rows = lv / 2;
        let cols = 1usize << (lv - log_rows);
        print_row(&measure(log_rows, cols, 32));
    }

    header("split: num_vars=16 fixed, vary log_rows, q=32 (Brakedown)");
    for &log_rows in &[2usize, 4, 6, 8, 10, 12] {
        let cols = 1usize << (16 - log_rows);
        print_row(&measure(log_rows, cols, 32));
    }

    header("queries: num_vars=16 square, vary t (Brakedown)");
    for &q in &[8usize, 16, 32, 64, 128, 256] {
        print_row(&measure(8, 1 << 8, q));
    }

    // Head-to-head: same shape & queries, only the linear code differs.
    header("Brakedown vs Ligero (RS rate 2): square split, q=32");
    for &lv in &[8usize, 12, 14, 16] {
        let log_rows = lv / 2;
        let cols = 1usize << (lv - log_rows);
        print_row(&measure_code("brakedown", bd_code(cols), log_rows, 32));
        print_row(&measure_code("ligero-rs", rs_code(cols), log_rows, 32));
    }
    println!(
        "\nRS distance δ=1−1/rate=0.5 vs Brakedown δ≈0.07 → at equal q each RS query is worth"
    );
    println!("~7x more security bits (RS needs ~7x fewer queries for the same soundness).");
    println!("Cost of that: RS codeword is denser AND RS encode is superlinear (see timing).");
}
