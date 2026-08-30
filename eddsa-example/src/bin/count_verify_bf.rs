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
        let c = p3_ed25519::counters::snapshot();
        // Local hash = Poseidon2-over-Ed25519 (counted via on_perm → PERM_COUNT); the BabyBear
        // config would use HASH_COUNT (`c.hash`) instead. Take whichever is populated.
        let perms = {
            let p = p3_ed25519::counters::PERM_COUNT.load(core::sync::atomic::Ordering::Relaxed);
            if p != 0 { p } else { c.hash }
        };
        let (mul, add) = (c.mul, c.add);
        // Default (ed25519 hash): its muls are in `mul`, so subtract them (perms × muls/perm) to
        // isolate the PCS ARITHMETIC muls. Under `bb-hash` the hash is BabyBear, so k=0 and
        // total mul == arithmetic mul (nothing to subtract).
        #[cfg(not(feature = "bb-hash"))]
        let k = eddsa_example::config::muls_per_perm();
        #[cfg(feature = "bb-hash")]
        let k = 0u64;
        let hash_muls = perms * k;
        let arith = mul.saturating_sub(hash_muls);
        println!(
            "PCS={:<16} num_bits={} ok={:<5} | total_muls={:>9} sq={:>9} inv={:>7} add={:>9} \
             perms={:>7} hash_muls={:>9} ARITH_muls={:>9}  (muls/perm={})",
            pcs, num_bits, res.is_ok(), mul, c.sq, c.inv, add, perms, hash_muls, arith, k
        );
    }
    #[cfg(not(feature = "count-ops"))]
    println!("PCS={pcs} verify_ok={:?} (build with --features count-ops to count field ops)", res.is_ok());
}
