//! M4: count ed25519 field ops in the eddsa multi-stark VERIFY, to validate Basefold's field-mul
//! advantage over the Brakedown tensor PCS. Build/run twice and compare:
//!   cargo run -p eddsa-example --release --features count-ops           --bin count_verify_bf
//!   cargo run -p eddsa-example --release --features count-ops,basefold  --bin count_verify_bf

use eddsa_example::config::{challenger, config_for_eddsa};
use eddsa_example::curve::{base_b, d, subgroup_order, valid_instance};
use p3_ed25519_air::eddsa::{
    EdDsaAir, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS, generate_trace_rows, public_values,
    trace_height,
};
use p3_multi_stark::{
    ProverInstance, ProverInstances, VerifierInstance, VerifierInstances, prove, setup, verify,
};
use p3_sumcheck::layout::Table;
use p3_util::log2_ceil_usize;

fn main() {
    let num_bits: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(253);
    let log_height = log2_ceil_usize(trace_height(num_bits));
    let (d_, base_b_, l) = (d(), base_b(), subgroup_order());
    let inst = valid_instance(num_bits, 3, 5, 7);
    let air = EdDsaAir { d: d_, base_b: base_b_, num_bits };
    let trace = generate_trace_rows(d_, base_b_, &l, inst.a, inst.r, &inst.s_bits, &inst.k_bits);
    let pis = public_values(inst.a, inst.r, &inst.s_bits, &inst.k_bits);
    let config = config_for_eddsa(log_height, NUM_EDDSA_COLS, NUM_EDDSA_PREPROCESSED_COLS);

    let (pk, vk) = setup(&config, &[&air], &mut challenger());
    let proof = prove(
        &config,
        ProverInstances::new(vec![ProverInstance::new(&air, Table::new(trace.transpose()), &pk, &pis)]),
        0,
        &mut challenger(),
    );

    #[cfg(feature = "count-ops")]
    p3_ed25519::counters::reset();

    let res = verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&air, &vk, log_height, &pis)]),
        &proof,
        0,
        &mut challenger(),
    );

    let pcs = if cfg!(feature = "basefold") { "Basefold" } else { "Brakedown/tensor" };
    #[cfg(feature = "count-ops")]
    {
        let (mul, add, perms) = p3_ed25519::counters::snapshot();
        // Locally the hash is Poseidon2-over-Ed25519, so its muls are in `mul`. Subtract them
        // (perms × muls/perm) to isolate the PCS ARITHMETIC muls — the metric the Basefold model
        // is about (on vultr2/guest the hash is BabyBear, so total mul == arithmetic mul).
        let k = eddsa_example::config::muls_per_perm();
        let hash_muls = perms * k;
        let arith = mul.saturating_sub(hash_muls);
        println!(
            "PCS={:<16} num_bits={} ok={:<5} | total_muls={:>9} perms={:>7} hash_muls={:>9} \
             ARITH_muls={:>9}  (muls/perm={})",
            pcs, num_bits, res.is_ok(), mul, perms, hash_muls, arith, k
        );
    }
    #[cfg(not(feature = "count-ops"))]
    println!("PCS={pcs} verify_ok={:?} (build with --features count-ops to count field ops)", res.is_ok());
}
