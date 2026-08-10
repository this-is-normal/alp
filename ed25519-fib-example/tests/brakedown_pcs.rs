//! End-to-end integration test: a native multilinear Fibonacci STARK over the
//! ed25519 base field, committed with the Brakedown tensor PCS.
//!
//! Establishes the two key properties:
//! 1. An honest Fibonacci proof round-trips through multi-stark + Brakedown over
//!    ed25519 (`F = EF = Ed25519`, a non-two-adic, non-`PrimeField64` field).
//! 2. A violated constraint is rejected by the zerocheck.
//!
//! The trace is kept small (`2^6` rows) because ed25519 arithmetic is slow.

use ed25519_fib_example::config::{challenger_brakedown, config_for_brakedown};
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
fn brakedown_pcs_honest_fibonacci_roundtrips_over_ed25519() {
    let n = 1 << LOG_HEIGHT;
    let trace = fib_trace(n);
    let pis = fib_public_values_for_trace(&trace);
    let config = config_for_brakedown(LOG_HEIGHT, NUM_COLS);

    let (pk, vk) = setup(&config, &[&FibAir], &mut challenger_brakedown());
    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(
            &FibAir,
            Table::new(trace.transpose()),
            &pk,
            &pis,
        )]),
        0,
        &mut challenger_brakedown(),
    );
    verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&FibAir, &vk, LOG_HEIGHT, &pis)]),
        &proof,
        0,
        &mut challenger_brakedown(),
    )
    .expect("honest Fibonacci proof must verify under the ed25519 Brakedown PCS");
}

#[test]
fn brakedown_pcs_rejects_violated_constraint_over_ed25519() {
    let n = 1 << LOG_HEIGHT;
    let mut trace = fib_trace(n);
    // Mutation: break a transition by shifting one row value the constraints read.
    trace.values[2 * NUM_COLS] += Ed25519::ONE;
    let pis = fib_public_values_for_trace(&trace);
    let config = config_for_brakedown(LOG_HEIGHT, NUM_COLS);

    let (pk, vk) = setup(&config, &[&FibAir], &mut challenger_brakedown());
    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(
            &FibAir,
            Table::new(trace.transpose()),
            &pk,
            &pis,
        )]),
        0,
        &mut challenger_brakedown(),
    );

    let err = verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&FibAir, &vk, LOG_HEIGHT, &pis)]),
        &proof,
        0,
        &mut challenger_brakedown(),
    )
    .expect_err("a broken constraint must be rejected");
    assert!(
        matches!(
            err,
            VerificationError::Zerocheck(ZerocheckError::FinalSumMismatch)
        ),
        "expected zerocheck final-sum mismatch, got {err:?}"
    );
}
