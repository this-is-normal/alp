//! Side-by-side end-to-end proof sizes for the same Fibonacci trace under each
//! PCS, so a regression in one is visible against the others.
//!
//! Not an assertion — a printed comparison:
//! `cargo test -p ed25519-fib-example --test pcs_comparison --release -- --nocapture`

use ed25519_fib_example::config::{
    challenger_brakedown, challenger_quasar, config_for_brakedown, config_for_quasar,
};
use ed25519_fib_example::fib_air::{FibAir, NUM_COLS, fib_public_values_for_trace, fib_trace};
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerifierInstance, VerifierInstances, prove, setup, verify,
};
use p3_stark_stats::serialized_size;
use p3_sumcheck::layout::Table;

const LOG_HEIGHT: usize = 6;

/// The generic bounds on `prove`/`verify` are long enough that a helper function
/// would be mostly `where` clauses; a macro keeps the comparison readable.
macro_rules! measure {
    ($name:literal, $config:expr, $challenger:expr) => {{
        let n = 1 << LOG_HEIGHT;
        let trace = fib_trace(n);
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
        let prove_us = t0.elapsed().as_micros();

        let t1 = std::time::Instant::now();
        verify(
            &config,
            VerifierInstances::new(vec![VerifierInstance::new(&FibAir, &vk, LOG_HEIGHT, &pis)]),
            &proof,
            0,
            &mut $challenger(),
        )
        .expect("honest proof must verify");
        let verify_us = t1.elapsed().as_micros();

        println!(
            "{:<10} proof {:>9} B   prove {:>8} us   verify {:>8} us",
            $name,
            serialized_size(&proof),
            prove_us,
            verify_us
        );
    }};
}

#[test]
fn proof_size_comparison() {
    println!("Fibonacci 2^{LOG_HEIGHT} rows x {NUM_COLS} cols over ed25519");
    measure!(
        "brakedown",
        config_for_brakedown(LOG_HEIGHT, NUM_COLS),
        challenger_brakedown
    );
    measure!(
        "quasar",
        config_for_quasar(LOG_HEIGHT, NUM_COLS),
        challenger_quasar
    );
}
