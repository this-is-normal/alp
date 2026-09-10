//! Generate a batched EdDSA proof and serialize `(num_bits, log_instances, pis,
//! vk, proof)` for the OpenVM recursion guest.
//!
//! ```text
//! cargo run -p eddsa-example --release --features "quasar,bb-hash,parallel" \
//!     --bin gen_proof -- 253 1 /root/eddsa_quasar.bin
//! ```
//!
//! The guest rebuilds the *same* config from `(log_height, width)` — every code
//! and permutation here is seeded deterministically — so only the proof, the
//! verifying key and the public values need to cross.

use eddsa_example::backend::NAME;
use eddsa_example::config::{challenger, config_for_eddsa};
use eddsa_example::curve::{base_b, d, valid_instances};
use eddsa_example::trace::generate_batch_trace;
use p3_ed25519_air::eddsa::{
    EdDsaAir, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS, public_values, trace_height,
};
use p3_multi_stark::{ProverInstance, ProverInstances, prove, setup};
use p3_sumcheck::layout::Table;
use p3_util::log2_ceil_usize;

fn main() {
    let num_bits: usize = std::env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(253);
    let num_signatures: usize = std::env::args()
        .nth(2)
        .and_then(|a| a.parse().ok())
        .unwrap_or(1);
    let path = std::env::args()
        .nth(3)
        .unwrap_or_else(|| "/root/eddsa_proof.bin".to_string());
    assert!(
        num_signatures.is_power_of_two(),
        "signature count must be a power of two, got {num_signatures}"
    );

    let log_instances = num_signatures.trailing_zeros() as usize;
    let instances = valid_instances(num_bits, num_signatures);
    let air = EdDsaAir { d: d(), base_b: base_b(), num_bits, log_instances };

    let trace = generate_batch_trace(&instances);
    let first = &instances[0];
    let pis = public_values(first.a, first.r, &first.s_bits, &first.k_bits);
    let log_height = log2_ceil_usize(trace_height(num_bits) * num_signatures);

    let config = config_for_eddsa(log_height, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS);
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

    // Quasar's coefficient commitments are public setup material: fixed by the
    // code seed, and the only thing the verifier reads out of a full
    // `QuasarPcs::new`. Shipping them (like the verifying key) lets the guest
    // call `QuasarPcs::new_verifier` and skip an O(N₂) encode-and-Merkle-commit
    // per proof. Every other backend ships nothing extra.
    #[cfg(feature = "quasar")]
    let quasar_setup = {
        let (main, preprocessed) = config.pcs_pair();
        // The resolved shape travels too: choosing it is an `f64` bisection over
        // transcendentals, which is free here and softfloat in a rv32 guest.
        (
            main.inner().coeff_commitment().clone(),
            preprocessed.inner().coeff_commitment().clone(),
            eddsa_example::backend::resolved_params(log_height, NUM_EDDSA_COLS),
            eddsa_example::backend::resolved_params(log_height, NUM_EDDSA_PREPROCESSED_COLS),
        )
    };
    #[cfg(not(feature = "quasar"))]
    let quasar_setup = ();

    let blob = postcard::to_allocvec(&(
        num_bits as u64,
        log_instances as u64,
        log_height as u64,
        &pis,
        &vk,
        &proof,
        &quasar_setup,
    ))
    .expect("serialize");
    std::fs::write(&path, &blob).expect("write");
    println!(
        "{NAME}: wrote {path}: {} bytes (num_bits={num_bits}, sigs={num_signatures}, \
         log_height={log_height})",
        blob.len()
    );
}
