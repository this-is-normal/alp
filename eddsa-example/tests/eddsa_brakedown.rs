//! End-to-end tests: the ed25519 `EdDsaAir` proved with the Brakedown multi-stark
//! prover over the native base field. Small `num_bits` keeps the trace tiny.

use eddsa_example::config::{F, challenger, config_for_eddsa_with_security};
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

/// This end-to-end test intentionally uses a tiny trace (fast). Its smallest committed poly is the
/// preprocessed one at ~2^10, where secure Brakedown at 100 bits is impossible (a 2^10 poly can open
/// at most ~1.7k columns, and 100-bit needs ≥2936). We therefore target a soundness level the trace
/// can honestly reach; `build_pcs` now *panics* rather than silently clamp to an insecure query count.
const TEST_SECURITY_BITS: usize = 40;

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
    let config = config_for_eddsa_with_security(
        log_height,
        NUM_EDDSA_COLS,
        NUM_EDDSA_PREPROCESSED_COLS,
        TEST_SECURITY_BITS,
    );
    (config, air, Table::new(trace.transpose()), pis, log_height)
}

#[test]
fn eddsa_brakedown_prove_and_verify() {
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
fn eddsa_brakedown_rejects_tampered_public_value() {
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
