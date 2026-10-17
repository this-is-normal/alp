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
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerifierInstance, VerifierInstances, prove, setup, verify,
};
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

    let config = config_for_eddsa(num_bits, log_height, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS);
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

    // `BLOB_STATS=1` breaks the blob down by component — the recursion guest
    // pays ~700 instructions per field element to parse this, so knowing which
    // part dominates says whether to attack the codec or the content.
    if std::env::var("BLOB_STATS").is_ok() {
        let n = |b: usize| format!("{b} B ({} field elements)", b / 32);
        eprintln!("  commitment          {}", n(p3_stark_stats::serialized_size(&proof.commitment)));
        eprintln!("  zerocheck sumcheck  {}", n(p3_stark_stats::serialized_size(&proof.sumcheck)));
        eprintln!("  main opening        {}", n(p3_stark_stats::serialized_size(&proof.opening)));
        eprintln!(
            "  preprocessed openin {}",
            n(p3_stark_stats::serialized_size(&proof.preprocessed_opening))
        );
        eprintln!("  verifying key       {}", n(p3_stark_stats::serialized_size(&vk)));
    }

    // Inverse hints for the recursion guest (strategy D).
    //
    // The guest needs `1/twist[idx]` wherever its queries land, but it cannot
    // know those positions until it runs — and we can, because verification is
    // deterministic. So: verify once here with the recorder on, harvest the
    // indices, and ship the inverses. The guest checks each with a single
    // multiply (`t · hint == 1`), so a wrong hint cannot be accepted.
    //
    // Without `record-queries` this is empty and the guest falls back to its
    // baked `TINV_WORDS` table, so the blob format does not fork.
    #[cfg(all(feature = "quasar", feature = "record-queries"))]
    let inverse_hints = {
        p3_basefold::query_recorder::reset();
        verify(
            &config,
            VerifierInstances::new(vec![VerifierInstance::new(&air, &vk, log_height, &pis)]),
            &proof,
            0,
            &mut challenger(),
        )
        .expect("honest proof must verify");
        let idxs = p3_basefold::query_recorder::take_sorted();
        let hints = eddsa_example::backend::twist_inverse_hints(&idxs);
        eprintln!(
            "  inverse hints: {} distinct indices, {} B",
            hints.0.len(),
            hints.0.len() * 36
        );
        hints
    };
    #[cfg(not(all(feature = "quasar", feature = "record-queries")))]
    let inverse_hints: (Vec<u32>, Vec<p3_ed25519::Ed25519>) = (Vec::new(), Vec::new());

    // Strategy M: commit to the twist table and ship the queried leaves with a
    // single pruned multi-path, instead of the guest carrying the whole table
    // in its program image (~86 MB of `.vmexe`, charged on every `prove app`).
    //
    // `MERKLE_TWISTS=1` turns it on. One tree covers both codes: they share a
    // seed and the index layout is a prefix, so a tree sized to the larger
    // `log_cols` contains every index the smaller can reach.
    #[cfg(all(feature = "quasar", feature = "record-queries"))]
    let merkle_twists = if std::env::var("MERKLE_TWISTS").is_ok() {
        let idxs = {
            p3_basefold::query_recorder::reset();
            verify(
                &config,
                VerifierInstances::new(vec![VerifierInstance::new(&air, &vk, log_height, &pis)]),
                &proof,
                0,
                &mut challenger(),
            )
            .expect("honest proof must verify");
            p3_basefold::query_recorder::take_sorted()
        };
        // One tree for both codes: same seed, prefix index layout.
        let k = eddsa_example::backend::resolved_params(log_height, NUM_EDDSA_COLS).0.log_cols.max(
            eddsa_example::backend::resolved_params(log_height, NUM_EDDSA_PREPROCESSED_COLS).0.log_cols,
        );
        let payload = eddsa_example::backend::twist_merkle_payload(k, &idxs);
        eprintln!(
            "  merkle twists: {} leaves, {} sibling digests ({} B), tree 2^{} leaves depth {}",
            payload.0.len(),
            payload.2.sibling_hashes.len(),
            payload.2.sibling_hashes.len() * 32,
            payload.3,
            payload.4,
        );
        Some((payload.0, payload.1, payload.2))
    } else {
        None
    };
    #[cfg(not(all(feature = "quasar", feature = "record-queries")))]
    let merkle_twists: Option<(
        Vec<u32>,
        Vec<p3_ed25519::Ed25519>,
        p3_merkle_tree::PrunedMerklePaths<p3_baby_bear::BabyBear, 8>,
    )> = None;

    // The flat format, not postcard: the guest deserialises this before it can
    // verify anything, and postcard's per-element machinery cost it ~710
    // instructions per field element against ~152 for a flat read. Slightly
    // larger on the wire (fixed-width lengths); the payload is incompressible
    // 32-byte field elements either way.
    let blob = p3_flat_serde::to_allocvec(&(
        num_bits as u64,
        log_instances as u64,
        log_height as u64,
        &pis,
        &vk,
        &proof,
        &quasar_setup,
        &inverse_hints,
        &merkle_twists,
    ))
    .expect("serialize");
    std::fs::write(&path, &blob).expect("write");
    println!(
        "{NAME}: wrote {path}: {} bytes (num_bits={num_bits}, sigs={num_signatures}, \
         log_height={log_height})",
        blob.len()
    );
}
