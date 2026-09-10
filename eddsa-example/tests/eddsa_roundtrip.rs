//! End-to-end ed25519 signature-verification STARK: prove and verify the real
//! `EdDsaAir` over the native ed25519 base field, and reject a tampered public
//! value.
//!
//! Backend-agnostic — runs against whichever commitment scheme the cargo
//! features select (Brakedown by default; `ligero`, `basefold` or `quasar`).


use eddsa_example::config::{F, challenger, config_for_eddsa};
use eddsa_example::curve::{base_b, d, subgroup_order, valid_instance};
use p3_ed25519_air::eddsa::{
    EdDsaAir, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS, generate_trace_rows, public_values,
    trace_height,
};
use p3_field::PrimeCharacteristicRing;
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerifierInstance, VerifierInstances, prove, setup, verify,
};
use p3_sumcheck::layout::Table;
use p3_util::log2_ceil_usize;

const NUM_BITS: usize = 32;

/// Build (config, air, trace, public values, log_height) for a valid instance.
fn fixture() -> (
    eddsa_example::config::EdDsaConfig,
    EdDsaAir<F>,
    Table<F>,
    Vec<F>,
    usize,
) {
    let log_height = log2_ceil_usize(trace_height(NUM_BITS));
    let inst = valid_instance(NUM_BITS, 3, 5, 7);
    let (d, base_b, l) = (d(), base_b(), subgroup_order());

    let air = EdDsaAir { d, base_b, num_bits: NUM_BITS, log_instances: 0 };
    let trace = generate_trace_rows(d, base_b, &l, inst.a, inst.r, &inst.s_bits, &inst.k_bits);
    let pis = public_values(inst.a, inst.r, &inst.s_bits, &inst.k_bits);
    let config = config_for_eddsa(log_height, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS);
    (config, air, Table::new(trace.transpose()), pis, log_height)
}

#[test]
fn eddsa_prove_and_verify() {
    let (config, air, table, pis, log_height) = fixture();
    let (pk, vk) = setup(&config, &[&air], &mut challenger());

    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(&air, table, &pk, &pis)]),
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
}

#[test]
fn eddsa_rejects_tampered_public_value() {
    let (config, air, table, pis, log_height) = fixture();
    let (pk, vk) = setup(&config, &[&air], &mut challenger());

    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(&air, table, &pk, &pis)]),
        0,
        &mut challenger(),
    );

    // Corrupt `s` (public value index 4): the verifier now checks a different claim
    // than the prover proved, so verification must fail.
    let mut bad_pis = pis;
    bad_pis[4] += F::ONE;

    let result = verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&air, &vk, log_height, &bad_pis)]),
        &proof,
        0,
        &mut challenger(),
    );
    assert!(result.is_err(), "a tampered public value must be rejected");
}

/// Batched proving: four signatures stacked into one table must round-trip, and
/// corrupting the *second* one must be caught — the property that makes stacking
/// meaningful. `EdDsaAir` gates its per-instance constraints on the preprocessed
/// `is_instance_first`/`is_instance_last` selectors rather than on the global
/// first row, so every stacked signature is fully constrained even though only
/// instance 0 is bound to the public values.
#[test]
fn eddsa_batched_roundtrips_and_rejects_a_corrupt_instance() {
    use eddsa_example::curve::valid_instances;
    use eddsa_example::trace::generate_batch_trace;
    use p3_ed25519_air::eddsa::{
        EdDsaAir, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS, public_values, trace_height,
    };

    const NUM_BITS: usize = 32;
    const NUM_SIGS: usize = 4;

    let instances = valid_instances(NUM_BITS, NUM_SIGS);
    let air = EdDsaAir {
        d: d(),
        base_b: base_b(),
        num_bits: NUM_BITS,
        log_instances: NUM_SIGS.trailing_zeros() as usize,
    };
    let first = &instances[0];
    let pis = public_values(first.a, first.r, &first.s_bits, &first.k_bits);
    let log_height = log2_ceil_usize(trace_height(NUM_BITS) * NUM_SIGS);
    let config = config_for_eddsa(log_height, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS);
    let (pk, vk) = setup(&config, &[&air], &mut challenger());

    let honest = generate_batch_trace(&instances);
    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(
            &air,
            Table::new(honest.transpose()),
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
    .expect("an honest batch of four signatures must verify");

    // Break a ladder row belonging to the SECOND instance.
    let mut bad = generate_batch_trace(&instances);
    let second_instance_row = trace_height(NUM_BITS) + 1;
    bad.values[second_instance_row * NUM_EDDSA_COLS] += F::ONE;
    let bad_proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(
            &air,
            Table::new(bad.transpose()),
            &pk,
            &pis,
        )]),
        0,
        &mut challenger(),
    );
    let result = verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&air, &vk, log_height, &pis)]),
        &bad_proof,
        0,
        &mut challenger(),
    );
    assert!(
        result.is_err(),
        "a corrupt second instance must be rejected — otherwise batching proves nothing about it"
    );
}
