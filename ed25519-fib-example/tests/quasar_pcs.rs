//! End-to-end integration test: a native multilinear Fibonacci STARK over the
//! ed25519 base field, committed with the **Quasar** PCS (QA code + code
//! switching into BaseFold).
//!
//! The Brakedown counterpart is `brakedown_pcs.rs`; this establishes the same
//! two properties through the polylog-verifier scheme.

use ed25519_fib_example::config::{challenger_quasar, config_for_quasar};
use ed25519_fib_example::fib_air::{FibAir, NUM_COLS, fib_public_values_for_trace, fib_trace};
use p3_ed25519::Ed25519;
use p3_field::PrimeCharacteristicRing;
use p3_multi_stark::zerocheck::ZerocheckError;
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerificationError, VerifierInstance, VerifierInstances, prove,
    setup, verify,
};
use p3_sumcheck::layout::Table;

const LOG_HEIGHT: usize = 6;

#[test]
fn quasar_pcs_honest_fibonacci_roundtrips_over_ed25519() {
    let n = 1 << LOG_HEIGHT;
    let trace = fib_trace(n);
    let pis = fib_public_values_for_trace(&trace);
    let config = config_for_quasar(LOG_HEIGHT, NUM_COLS);

    let (pk, vk) = setup(&config, &[&FibAir], &mut challenger_quasar());
    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(
            &FibAir,
            Table::new(trace.transpose()),
            &pk,
            &pis,
        )]),
        0,
        &mut challenger_quasar(),
    );
    verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&FibAir, &vk, LOG_HEIGHT, &pis)]),
        &proof,
        0,
        &mut challenger_quasar(),
    )
    .expect("honest Fibonacci proof must verify under the ed25519 Quasar PCS");
}

#[test]
fn quasar_pcs_rejects_violated_constraint_over_ed25519() {
    let n = 1 << LOG_HEIGHT;
    let mut trace = fib_trace(n);
    trace.values[2 * NUM_COLS] += Ed25519::ONE;
    let pis = fib_public_values_for_trace(&trace);
    let config = config_for_quasar(LOG_HEIGHT, NUM_COLS);

    let (pk, vk) = setup(&config, &[&FibAir], &mut challenger_quasar());
    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(
            &FibAir,
            Table::new(trace.transpose()),
            &pk,
            &pis,
        )]),
        0,
        &mut challenger_quasar(),
    );
    let err = verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&FibAir, &vk, LOG_HEIGHT, &pis)]),
        &proof,
        0,
        &mut challenger_quasar(),
    )
    .expect_err("a violated transition must be rejected");
    assert!(
        matches!(
            err,
            VerificationError::Zerocheck(ZerocheckError::FinalSumMismatch)
        ),
        "expected zerocheck final-sum mismatch, got {err:?}"
    );
}
