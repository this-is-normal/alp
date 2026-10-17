//! Integration tests for the trivial field-agnostic PCS on the multi-stark path.
//!
//! These establish the two properties M1 is meant to demonstrate:
//! 1. An honest multilinear Fibonacci proof round-trips through multi-stark
//!    with a PCS that requires no two-adicity (field-agnostic proving works).
//! 2. multi-stark's zerocheck still rejects a constraint violation even with
//!    the trivial PCS — i.e. the trivial PCS integrates correctly, it does not
//!    accidentally accept everything.

use multi_fib_example::config::{challenger_trivial, config_for_trivial};
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
fn trivial_pcs_honest_fibonacci_roundtrips() {
    let n = 1 << LOG_HEIGHT;
    let trace = fib_trace(n);
    let pis = fib_public_values_for_trace(&trace);
    let config = config_for_trivial(LOG_HEIGHT, NUM_COLS);

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
        VerifierInstances::new(vec![VerifierInstance::new(&FibAir, &vk, LOG_HEIGHT, &pis)]),
        &proof,
        0,
        &mut challenger_trivial(),
    )
    .expect("honest Fibonacci proof must verify under the trivial PCS");
}

#[test]
fn trivial_pcs_rejects_violated_constraint() {
    let n = 1 << LOG_HEIGHT;
    let mut trace = fib_trace(n);
    // Mutation: break a transition by shifting one row value the constraints read.
    trace.values[2 * NUM_COLS] += BabyBear::ONE;
    let pis = fib_public_values_for_trace(&trace);
    let config = config_for_trivial(LOG_HEIGHT, NUM_COLS);

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

    // The zerocheck closes on a nonzero constraint value regardless of the PCS.
    let err = verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&FibAir, &vk, LOG_HEIGHT, &pis)]),
        &proof,
        0,
        &mut challenger_trivial(),
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
