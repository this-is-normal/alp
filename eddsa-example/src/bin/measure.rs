//! Measure the end-to-end EdDSA proof under whichever commitment backend the
//! cargo features select: proof size and prove/verify wall-clock.
//!
//! ```text
//! for f in "" ligero basefold quasar; do
//!   cargo run -p eddsa-example --release --bin measure ${f:+--features $f} -- 253 1,4,16,64,256
//! done
//! ```
//! Arg 1 is the scalar bit length (default 253, the real ed25519 size). Arg 2 is a
//! comma-separated list of signature counts to batch into one proof (powers of two,
//! default `1`) — `EdDsaAir` stacks them vertically, so the verifier does the work of
//! one signature plus `log₂n` extra sumcheck rounds.

use eddsa_example::backend::NAME;
use eddsa_example::config::{challenger, config_for_eddsa};
use eddsa_example::curve::valid_instances;
use eddsa_example::trace::generate_batch_trace;
use p3_ed25519_air::eddsa::{
    EdDsaAir, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS, public_values, trace_height,
};
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerifierInstance, VerifierInstances, prove, setup, verify,
};
use p3_stark_stats::serialized_size;
use p3_sumcheck::layout::Table;
use p3_util::log2_ceil_usize;
use tracing_forest::ForestLayer;
use tracing_forest::util::LevelFilter;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, Registry};

fn main() {
    // `PROFILE=1` prints the prover's span tree (per-phase wall-clock) instead
    // of just the summary line — the only way to see where a phase-level
    // parallel speedup is being lost.
    if std::env::var("PROFILE").is_ok() {
        Registry::default()
            .with(
                EnvFilter::builder()
                    .with_default_directive(LevelFilter::INFO.into())
                    .from_env_lossy(),
            )
            .with(ForestLayer::default())
            .init();
    }

    let num_bits: usize = std::env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(253);
    let counts: Vec<usize> = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "1".to_string())
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    for num_signatures in counts {
        assert!(
            num_signatures.is_power_of_two(),
            "signature count must be a power of two, got {num_signatures}"
        );
        measure(num_bits, num_signatures);
    }
}

fn measure(num_bits: usize, num_signatures: usize) {
    let log_instances = num_signatures.trailing_zeros() as usize;
    let instances = valid_instances(num_bits, num_signatures);
    let (d, base_b) = (eddsa_example::curve::d(), eddsa_example::curve::base_b());
    let air = EdDsaAir { d, base_b, num_bits, log_instances };

    let trace = generate_batch_trace(&instances);
    // Only instance 0 is bound to the public values (see `trace`'s scope note).
    let first = &instances[0];
    let pis = public_values(first.a, first.r, &first.s_bits, &first.k_bits);
    let height = trace_height(num_bits) * num_signatures;
    let log_height = log2_ceil_usize(height);

    let t_setup = std::time::Instant::now();
    let config = config_for_eddsa(num_bits, log_height, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS);
    let (pk, vk) = setup(&config, &[&air], &mut challenger());
    let setup_ms = t_setup.elapsed().as_secs_f64() * 1e3;

    let t_prove = std::time::Instant::now();
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
    let prove_ms = t_prove.elapsed().as_secs_f64() * 1e3;

    let t_verify = std::time::Instant::now();
    verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&air, &vk, log_height, &pis)]),
        &proof,
        0,
        &mut challenger(),
    )
    .expect("honest EdDSA proof must verify");
    let verify_ms = t_verify.elapsed().as_secs_f64() * 1e3;

    let bytes = serialized_size(&proof);
    println!(
        "{NAME:<10} sigs {num_signatures:>4}  trace={height:>6}x{NUM_EDDSA_COLS}  \
         proof {:>9} B ({:>7.1} B/sig)  setup {setup_ms:>8.1} ms  \
         prove {prove_ms:>9.1} ms ({:>8.2} ms/sig)  verify {verify_ms:>8.1} ms",
        bytes,
        bytes as f64 / num_signatures as f64,
        prove_ms / num_signatures as f64,
    );
}
