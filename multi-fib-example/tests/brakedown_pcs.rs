//! Integration tests for the Brakedown/Ligero tensor PCS on the multi-stark path.
//!
//! These mirror the trivial-PCS tests and establish the two M3b properties:
//! 1. An honest multilinear Fibonacci proof round-trips through multi-stark when
//!    committed with the Brakedown tensor PCS (the stacked-layout sumcheck
//!    reduction plus the tensor point opening replace WHIR end to end).
//! 2. multi-stark's zerocheck still rejects a constraint violation under the
//!    Brakedown PCS — the PCS integrates correctly and does not accept
//!    everything.

use multi_fib_example::config::{challenger_brakedown, config_for_brakedown};
use multi_fib_example::fib_air::{FibAir, NUM_COLS, fib_public_values_for_trace, fib_trace};
use p3_baby_bear::BabyBear;
use p3_field::PrimeCharacteristicRing;
use p3_multi_stark::zerocheck::ZerocheckError;
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerificationError, VerifierInstance, VerifierInstances, prove,
    setup, verify,
};
use p3_sumcheck::layout::Table;

const LOG_HEIGHT: usize = 8;

#[test]
fn brakedown_pcs_honest_fibonacci_roundtrips() {
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
    .expect("honest Fibonacci proof must verify under the Brakedown PCS");
}

#[test]
fn brakedown_pcs_rejects_violated_constraint() {
    let n = 1 << LOG_HEIGHT;
    let mut trace = fib_trace(n);
    // Mutation: break a transition by shifting one row value the constraints read.
    trace.values[2 * NUM_COLS] += BabyBear::ONE;
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

    // The zerocheck closes on a nonzero constraint value regardless of the PCS.
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
