//! How many foldable-code positions does a verifier actually read?
//!
//! This is the number that sizes both candidate replacements for the baked
//! twist table:
//!
//! - **hinted inverses** ship one field element per *distinct* index, so
//!   `distinct x 32 B` is the blob cost;
//! - **a Merkle-committed table** ships a leaf *and a path* per distinct index,
//!   so the blob cost is `distinct x (32 + depth x 32) B`.
//!
//! and `reads` (with repeats) is what predicts guest instruction cost, since
//! the guest pays per read, not per distinct element.
//!
//! ```text
//! cargo run -p eddsa-example --release --features "quasar,bb-hash,parallel,record-queries" \
//!     --bin count_queries -- 253 32,64,128,256,512,1024
//! ```

use eddsa_example::backend::NAME;
use eddsa_example::config::{challenger, config_for_eddsa};
use eddsa_example::curve::valid_instances;
use eddsa_example::trace::generate_batch_trace;
use p3_basefold::query_recorder;
use p3_ed25519_air::eddsa::{
    EdDsaAir, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS, public_values, trace_height,
};
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerifierInstance, VerifierInstances, prove, setup, verify,
};
use p3_sumcheck::layout::Table;
use p3_util::log2_ceil_usize;

/// Must match the backend's own foldable-code rate.
const BF_INV_RATE: usize = 4;

fn main() {
    let num_bits: usize = std::env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(253);
    let counts: Vec<usize> = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "256".to_string())
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    println!(
        "{:<8} {:>5} {:>9} {:>9} {:>8} {:>7} {:>12} {:>12} {:>10}",
        "backend",
        "sigs",
        "reads",
        "distinct",
        "maxidx",
        "levels",
        "hint_bytes",
        "merkle_bytes",
        "depth"
    );
    for n in counts {
        assert!(n.is_power_of_two(), "signature count must be a power of two");
        run(num_bits, n);
    }
}

fn run(num_bits: usize, num_signatures: usize) {
    let log_instances = num_signatures.trailing_zeros() as usize;
    let instances = valid_instances(num_bits, num_signatures);
    let air = EdDsaAir {
        d: eddsa_example::curve::d(),
        base_b: eddsa_example::curve::base_b(),
        num_bits,
        log_instances,
    };
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

    // Only the verify is recorded — the prover touches every position, which
    // says nothing about what a recursion guest needs.
    query_recorder::reset();
    verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&air, &vk, log_height, &pis)]),
        &proof,
        0,
        &mut challenger(),
    )
    .expect("honest EdDSA proof must verify");

    let reads = query_recorder::read_count();
    let idxs = query_recorder::take_sorted();
    let distinct = idxs.len();
    let max_idx = idxs.last().copied().unwrap_or(0);

    // Highest level any index lands in: level L spans [c·2^{L-1}, c·2^L).
    let levels = (1..=32)
        .find(|&l| max_idx < BF_INV_RATE << l)
        .unwrap_or(0);
    // A tree over the whole index space the verifier can reach.
    let leaves = (BF_INV_RATE << levels).next_power_of_two();
    let depth = log2_ceil_usize(leaves);

    // 32 B per ed25519 element; a BabyBear digest is 8 lanes = 32 B.
    let hint_bytes = distinct * 32;
    let merkle_bytes = distinct * (32 + depth * 32);

    println!(
        "{NAME:<8} {num_signatures:>5} {reads:>9} {distinct:>9} {max_idx:>8} {levels:>7} \
         {hint_bytes:>12} {merkle_bytes:>12} {depth:>10}"
    );
}
