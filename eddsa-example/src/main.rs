//! Prove and verify ed25519 signature verification (`EdDsaAir`) over the native
//! ed25519 base field with the Brakedown multi-stark prover.
//!
//! Run (from the workspace root, into a LOCAL target dir — sshfs can't link):
//! ```text
//! CARGO_TARGET_DIR=/home/ubuntu/alp-target cargo run -p eddsa-example --release
//! CARGO_TARGET_DIR=/home/ubuntu/alp-target cargo run -p eddsa-example --release -- 64
//! ```
//! Arg 1 is the scalar bit length `num_bits` (default 253, the real ed25519 size);
//! the trace is `next_power_of_two(num_bits + 3)` rows tall.

use eddsa_example::config::{F, challenger, config_for_eddsa};
use eddsa_example::curve::{base_b, d, subgroup_order, valid_instance};
use p3_ed25519_air::eddsa::{
    EdDsaAir, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS, generate_trace_rows, public_values,
    trace_height,
};
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerifierInstance, VerifierInstances, prove, setup, verify,
};
use p3_sumcheck::layout::Table;
use p3_util::log2_ceil_usize;
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

    let num_bits: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(253);

    let height = trace_height(num_bits);
    let log_height = log2_ceil_usize(height);

    tracing::info!(
        "EdDSA verification over ed25519: num_bits = {num_bits}, trace = {height} rows x {NUM_EDDSA_COLS} cols \
         (+ {NUM_EDDSA_PREPROCESSED_COLS} preprocessed cols)"
    );

    // A genuine verifying instance (small secrets a=3, r=5, k=7 ⇒ s = 26); the
    // circuit is full `num_bits` rows regardless of scalar magnitude. Constructing
    // this also independently checks R + [k]A == [s]B over the real curve.
    let instance = valid_instance(num_bits, 3, 5, 7);
    let (d, base_b, l) = (d(), base_b(), subgroup_order());

    let air = EdDsaAir { d, base_b, num_bits };
    let trace = generate_trace_rows(
        d,
        base_b,
        &l,
        instance.a,
        instance.r,
        &instance.s_bits,
        &instance.k_bits,
    );
    let pis: Vec<F> = public_values(instance.a, instance.r, &instance.s_bits, &instance.k_bits);

    let config = config_for_eddsa(log_height, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS);

    // Setup commits the (witness-independent) preprocessed trace once.
    let (pk, vk) = setup(&config, &[&air], &mut challenger());

    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(
            &air,
            Table::new(trace.transpose()),
            &pk,
            &pis,
        )]),
        0,
        &mut challenger(),
    );

    verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&air, &vk, log_height, &pis)]),
        &proof,
        0,
        &mut challenger(),
    )
    .expect("honest EdDSA proof must verify");

    println!(
        "OK: ed25519 signature verification proved and verified over the native base field \
         (num_bits = {num_bits}, {height} rows)."
    );
}
