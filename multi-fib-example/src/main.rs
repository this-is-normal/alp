//! Minimal runnable **multilinear** Fibonacci proof: multi-stark (SuperSpartan-
//! style zerocheck) with a pluggable multilinear PCS, over BabyBear.
//!
//! This is the multivariate-path counterpart to `fib-example` (the univariate
//! uni-stark/FRI path), and the M0/M1 harness for the Brakedown work:
//!
//! * `whir`    (default) — commit with the WHIR PCS (needs two-adicity).
//! * `trivial`           — commit with the field-agnostic trivial oracle PCS
//!                         (no two-adicity), proving the path is field-agnostic.
//!
//! Run (from the workspace root, into a LOCAL target dir — sshfs can't link):
//!
//! ```text
//! CARGO_TARGET_DIR=/home/ubuntu/alp-target cargo run -p multi-fib-example --release -- 8
//! CARGO_TARGET_DIR=/home/ubuntu/alp-target cargo run -p multi-fib-example --release -- 8 trivial
//! ```
//!
//! Arg 1 is `log2(rows)` (default 8 → 256 rows); arg 2 is the PCS (default `whir`).

use multi_fib_example::config::{
    F, MultiFibConfig, MultiFibConfigTrivial, challenger, challenger_trivial, config_for,
    config_for_trivial,
};
use multi_fib_example::fib_air::{FibAir, NUM_COLS, fib_public_values_for_trace, fib_trace};
use p3_matrix::dense::RowMajorMatrix;
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerifierInstance, VerifierInstances, prove, setup, verify,
};
use p3_sumcheck::layout::Table;
use tracing_forest::ForestLayer;
use tracing_forest::util::LevelFilter;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, Registry};

fn main() {
    let env_filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .from_env_lossy();
    Registry::default()
        .with(env_filter)
        .with(ForestLayer::default())
        .init();

    let log_height: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);
    let pcs_kind = std::env::args().nth(2).unwrap_or_else(|| "whir".into());
    let n = 1 << log_height;

    let trace = fib_trace(n);
    let pis = fib_public_values_for_trace(&trace);

    tracing::info!(
        "multilinear Fibonacci [{pcs_kind}]: rows = 2^{log_height} = {n}, fib({}) = {}",
        n - 1,
        pis[2]
    );

    match pcs_kind.as_str() {
        "trivial" => run_trivial(log_height, trace, pis),
        "whir" => run_whir(log_height, trace, pis),
        other => panic!("unknown pcs '{other}': expected 'whir' or 'trivial'"),
    }

    println!("OK: multilinear Fibonacci ({pcs_kind}) proved and verified over 2^{log_height} = {n} rows.");
}

/// Prove + verify with the WHIR PCS.
fn run_whir(log_height: usize, trace: RowMajorMatrix<F>, pis: [F; 3]) {
    let config: MultiFibConfig = config_for(log_height, NUM_COLS);

    let (pk, vk) = setup(&config, &[&FibAir], &mut challenger(&config));
    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(
            &FibAir,
            Table::new(trace.transpose()),
            &pk,
            &pis,
        )]),
        0,
        &mut challenger(&config),
    );
    verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&FibAir, &vk, log_height, &pis)]),
        &proof,
        0,
        &mut challenger(&config),
    )
    .expect("honest Fibonacci proof must verify (whir)");
}

/// Prove + verify with the trivial field-agnostic PCS.
fn run_trivial(log_height: usize, trace: RowMajorMatrix<F>, pis: [F; 3]) {
    let config: MultiFibConfigTrivial = config_for_trivial(log_height, NUM_COLS);

    let (pk, vk) = setup(&config, &[&FibAir], &mut challenger_trivial());
    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(
            &FibAir,
            Table::new(trace.transpose()),
            &pk,
            &pis,
        )]),
        0,
        &mut challenger_trivial(),
    );
    verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&FibAir, &vk, log_height, &pis)]),
        &proof,
        0,
        &mut challenger_trivial(),
    )
    .expect("honest Fibonacci proof must verify (trivial)");
}
