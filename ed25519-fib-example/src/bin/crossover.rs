//! Where does Quasar's polylogarithmic proof overtake Brakedown's square-root
//! one? Sweep the trace height and print both, end to end.
//!
//! ```text
//! cargo run -p ed25519-fib-example --release --bin crossover -- 6,10,14,16
//! ```
//!
//! Quasar carries a large *fixed* code-switching cost (two BaseFold openings),
//! so it loses at toy sizes and gains as `t·N₁` — Brakedown's square-root term —
//! grows.
//!
//! For the authoritative answer use `eddsa-example --bin measure`, which sweeps
//! signature counts on the real AIR with both backends at their real query
//! counts; this Fibonacci sweep is a quick sanity check on a 2-column trace.

use ed25519_fib_example::config::{
    challenger_brakedown, challenger_quasar, config_for_brakedown, config_for_quasar,
};
use ed25519_fib_example::fib_air::{FibAir, NUM_COLS, fib_public_values_for_trace, fib_trace};
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerifierInstance, VerifierInstances, prove, setup, verify,
};
use p3_stark_stats::serialized_size;
use p3_sumcheck::layout::Table;

/// The generic bounds on `prove`/`verify` make a helper mostly `where` clauses.
macro_rules! measure {
    ($name:literal, $log_height:expr, $config:expr, $challenger:expr) => {{
        let log_height = $log_height;
        let trace = fib_trace(1 << log_height);
        let pis = fib_public_values_for_trace(&trace);
        let config = $config;

        let (pk, vk) = setup(&config, &[&FibAir], &mut $challenger());
        let t0 = std::time::Instant::now();
        let proof = prove(
            &config,
            ProverInstances::new(vec![ProverInstance::new(
                &FibAir,
                Table::new(trace.transpose()),
                &pk,
                &pis,
            )]),
            0,
            &mut $challenger(),
        );
        let prove_ms = t0.elapsed().as_secs_f64() * 1e3;

        let t1 = std::time::Instant::now();
        verify(
            &config,
            VerifierInstances::new(vec![VerifierInstance::new(&FibAir, &vk, log_height, &pis)]),
            &proof,
            0,
            &mut $challenger(),
        )
        .expect("honest proof must verify");
        let verify_ms = t1.elapsed().as_secs_f64() * 1e3;

        println!(
            "  {:<10} proof {:>9} B   prove {:>9.1} ms   verify {:>8.1} ms",
            $name,
            serialized_size(&proof),
            prove_ms,
            verify_ms
        );
        serialized_size(&proof)
    }};
}

fn main() {
    let heights: Vec<usize> = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "6,10,12,14".to_string())
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    for log_height in heights {
        println!("trace 2^{log_height} x {NUM_COLS} over ed25519 (stacked arity {})", log_height + 1);
        let bd = measure!(
            "brakedown",
            log_height,
            config_for_brakedown(log_height, NUM_COLS),
            challenger_brakedown
        );
        let qz = measure!(
            "quasar",
            log_height,
            config_for_quasar(log_height, NUM_COLS),
            challenger_quasar
        );
        println!(
            "  → quasar/brakedown proof ratio {:.2}x{}",
            qz as f64 / bd as f64,
            if qz < bd { "  (quasar ahead)" } else { "" }
        );
    }
}
