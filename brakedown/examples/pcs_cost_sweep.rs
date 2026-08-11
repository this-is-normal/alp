//! Brakedown tensor-PCS **cost-model sweep** (deterministic, not timed).
//!
//! Drives the real recursive Brakedown `TensorPcs` through the *instrumented*
//! hash stack (`p3-stark-stats` counters on leaf hashing, node compression, and
//! the challenger sponge) and prints, per configuration, the exact prover- and
//! verifier-side hash-op counts plus the serialized proof size. These are the
//! inputs to the cost model
//!
//! ```text
//! time ≈ Σ_role  op_count[role] · per_op_time[role]
//! ```
//!
//! where `op_count` comes from *this* sweep (deterministic, machine-independent)
//! and `per_op_time` comes from the criterion microbenches. Counting is kept out
//! of any timed path — the plain `TensorPcs` (see `benches/tensor_pcs.rs`) is
//! what gets timed; this binary only counts.
//!
//! Three sub-sweeps expose the three tuning knobs:
//!   1. **scale**  — grow `num_vars` (square-ish row/col split);
//!   2. **split**  — fix `num_vars`, slide the row/column boundary `log_rows`;
//!   3. **queries** — fix shape, vary the column-query count `t`.
//!
//! Run (into a LOCAL target dir — sshfs can't link):
//! ```text
//! CARGO_TARGET_DIR=/home/ubuntu/alp-target \
//!   cargo run -p p3-brakedown --example pcs_cost_sweep --release
//! ```
//!
//! NOTE: BabyBear uses *packed* (SIMD) Merkle hashing, so one counted leaf-hash
//! invocation covers `F::Packing::WIDTH` leaves. The counts are hash
//! *invocations*, not per-leaf; the `compress`/`sponge` columns are unaffected.
//! A `Packing = Self` field (bn254/ed25519) would make the leaf column exact.

use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_brakedown::{BrakedownCode, LinearCode, TensorPcs};
use p3_challenger::DuplexChallenger;
use p3_field::extension::BinomialExtensionField;
use p3_field::Field;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_multilinear_util::point::Point;
use p3_stark_stats::{serialized_size, CountingCompress, CountingHasher, CountingPerm, HashCounts};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

type F = BabyBear;
type EF = BinomialExtensionField<F, 4>;
type Perm = Poseidon2BabyBear<16>;
type MyHash = CountingHasher<PaddingFreeSponge<Perm, 16, 8, 8>>;
type MyCompress = CountingCompress<TruncatedPermutation<Perm, 2, 8, 16>>;
type PackedF = <F as Field>::Packing;
type MyMmcs = MerkleTreeMmcs<PackedF, PackedF, MyHash, MyCompress, 2, 8>;
type MyChallenger = DuplexChallenger<F, CountingPerm<Perm>, 16, 8>;
type Pcs = TensorPcs<F, EF, BrakedownCode<F, EF>, MyMmcs>;

/// Fixed seeds so every row is reproducible and the two sides share one code.
const CODE_SEED: u64 = 0xBADC0DE;
const PERM_SEED: u64 = 0xC0FFEE;

fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(PERM_SEED);
    Perm::new_from_rng_128(&mut rng)
}

/// One measured configuration.
struct Row {
    num_vars: usize,
    rows: usize,
    cols: usize,
    code_n: usize,
    queries: usize,
    prover: (usize, usize, usize),
    verifier: (usize, usize, usize),
    proof_bytes: usize,
}

/// Build an instrumented Brakedown tensor PCS for `(log_rows, cols)`, sharing the
/// three op counters so prover/verifier work can be attributed separately.
fn build_pcs(log_rows: usize, cols: usize, num_queries: usize, counts: &HashCounts) -> Pcs {
    let hash = CountingHasher::with_counter(PaddingFreeSponge::new(perm()), counts.leaf_hash.clone());
    let compress =
        CountingCompress::with_counter(TruncatedPermutation::new(perm()), counts.compress.clone());
    let mmcs = MyMmcs::new(hash, compress, 0);

    let mut code_rng = SmallRng::seed_from_u64(CODE_SEED);
    let code = BrakedownCode::<F, EF>::from_message_len(cols, &mut code_rng);
    debug_assert!(num_queries <= code.codeword_len(), "queries exceed codeword length");
    TensorPcs::new(code, mmcs, log_rows, num_queries)
}

fn counting_challenger(counts: &HashCounts) -> MyChallenger {
    DuplexChallenger::new(CountingPerm::with_counter(perm(), counts.sponge.clone()))
}

/// Run one full commit→open→verify roundtrip under the counters, attributing
/// commit+open to the prover and verify to the verifier.
fn measure(log_rows: usize, cols: usize, num_queries: usize) -> Row {
    let num_vars = log_rows + cols.trailing_zeros() as usize;
    let counts = HashCounts::new();
    let code_n = {
        // Read the codeword length first, so we can clamp the query count to it.
        let mut r = SmallRng::seed_from_u64(CODE_SEED);
        BrakedownCode::<F, EF>::from_message_len(cols, &mut r).codeword_len()
    };
    // A proof can query at most every column; clamp (tiny shapes have n < t).
    let num_queries = num_queries.min(code_n);
    let pcs = build_pcs(log_rows, cols, num_queries, &counts);

    let mut rng = SmallRng::seed_from_u64(1);
    let evals: Vec<F> = (0..(1usize << num_vars)).map(|_| rng.random::<F>()).collect();
    let point = Point::new((0..num_vars).map(|_| rng.random::<EF>()).collect());

    // Prover: commit + open on one challenger.
    counts.reset();
    let mut ch_p = counting_challenger(&counts);
    let (commitment, prover_data) = pcs.commit(evals, &mut ch_p);
    let (value, proof) = pcs.open(&prover_data, &point, &mut ch_p);
    let prover = counts.snapshot();
    let proof_bytes = serialized_size(&proof);

    // Verifier: fresh identically-seeded challenger.
    counts.reset();
    let mut ch_v = counting_challenger(&counts);
    pcs.verify(&commitment, &point, value, &proof, &mut ch_v)
        .expect("honest opening must verify");
    let verifier = counts.snapshot();

    Row {
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

fn header(title: &str) {
    println!("\n=== {title} ===");
    println!(
        "{:>4} {:>6} {:>6} {:>7} {:>4} | {:>16} | {:>16} | {:>9}",
        "vars", "rows", "cols", "code_n", "q", "prover L/C/S", "verifier L/C/S", "proof(B)"
    );
}

fn print_row(r: &Row) {
    let (pl, pc, ps) = r.prover;
    let (vl, vc, vs) = r.verifier;
    println!(
        "{:>4} {:>6} {:>6} {:>7} {:>4} | {:>16} | {:>16} | {:>9}",
        r.num_vars,
        r.rows,
        r.cols,
        r.code_n,
        r.queries,
        format!("{pl}/{pc}/{ps}"),
        format!("{vl}/{vc}/{vs}"),
        r.proof_bytes,
    );
}

fn main() {
    println!("Brakedown tensor-PCS cost-model sweep (BabyBear + Poseidon2 + Merkle)");
    println!("counts are hash INVOCATIONS: L=leaf-hash  C=compress  S=sponge-perm");
    println!("(BabyBear packs leaves — L covers Packing::WIDTH leaves per call)");

    // 1) SCALE: grow num_vars, square-ish split, fixed query budget.
    header("scale: num_vars up, square split, q=32");
    for &lv in &[8usize, 12, 14, 16, 18, 20] {
        let log_rows = lv / 2;
        let cols = 1usize << (lv - log_rows);
        print_row(&measure(log_rows, cols, 32));
    }

    // 2) SPLIT: fix num_vars = 16, slide the row/column boundary. Wider rows
    //    (bigger cols) mean more encoding work but fewer/longer combined rows;
    //    this knob trades prover encode cost against proof size and verify cost.
    header("split: num_vars=16 fixed, vary log_rows, q=32");
    for &log_rows in &[2usize, 4, 6, 8, 10, 12] {
        let cols = 1usize << (16 - log_rows);
        print_row(&measure(log_rows, cols, 32));
    }

    // 3) QUERIES: fix a square 2^16 shape, vary the column-query count t. This
    //    is the soundness knob: t ≈ ⌈−λ / log2(1−δ/3)⌉ sets the security level.
    header("queries: num_vars=16 square, vary t");
    for &q in &[8usize, 16, 32, 64, 128, 256] {
        print_row(&measure(8, 1 << 8, q));
    }

    println!("\nprover  cost ~ commit(encode+Merkle build) + open(t Merkle paths)");
    println!("verifier cost ~ t Merkle-path checks (≈ t · tree_depth compressions)");
}
