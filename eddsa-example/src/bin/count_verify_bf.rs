//! Attribute the VERIFY cost between **hashing** and **field arithmetic**.
//!
//! Counts, for one batched EdDSA verify: ed25519 field ops, Poseidon2
//! permutations, and cross-field limb splits (the glue every hashed element
//! goes through). Multiply each by its *measured in-guest* instruction cost
//! (see `alp-openvm/guests/microbench`) to get where an OpenVM recursion
//! verifier actually spends itself.
//!
//! ```text
//! for f in "quasar" "basefold" ""; do
//!   cargo run -p eddsa-example --release --features "bb-hash,count-ops${f:+,$f}" \
//!       --bin count_verify_bf -- 253 16
//! done
//! ```

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

/// Measured cost of one operation inside the OpenVM guest, in RISC-V
/// instructions (`alp-openvm/guests/microbench`, rv32im + modular + poseidon2
/// precompiles). These are what turn op counts into a time attribution.
#[cfg(feature = "count-ops")]
mod guest_cost {
    pub const PERM: u64 = 16;
    pub const MUL: u64 = 186;
    pub const SQ: u64 = 186;
    pub const ADD: u64 = 55;
    /// Measured, not modelled: 14,879, well under the 252·SQ + 77·MUL a Fermat
    /// inversion would cost, because the guest routes it through OpenVM's
    /// modular precompile rather than exponentiating.
    pub const INV: u64 = 14_879;
    pub const SPLIT: u64 = 222;
}

fn main() {
    let num_bits: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(253);
    let num_signatures: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(1);
    assert!(num_signatures.is_power_of_two());

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

    #[cfg(feature = "count-ops")]
    {
        p3_ed25519::counters::reset();
        p3_field::split_counter::reset();
    }

    let res = verify(
        &config,
        VerifierInstances::new(vec![VerifierInstance::new(&air, &vk, log_height, &pis)]),
        &proof,
        0,
        &mut challenger(),
    );

    #[cfg(feature = "count-ops")]
    {
        use core::sync::atomic::Ordering;
        let c = p3_ed25519::counters::snapshot();
        let perms = {
            let p = p3_ed25519::counters::PERM_COUNT.load(Ordering::Relaxed);
            if p != 0 { p } else { c.hash }
        };
        let splits = p3_field::split_counter::snapshot();

        // Hashing = the permutations plus the limb decomposition that feeds them.
        let perm_cost = perms * guest_cost::PERM;
        let split_cost = splits * guest_cost::SPLIT;
        let hash_cost = perm_cost + split_cost;
        // Arithmetic = ed25519 field ops. An inversion is a^(p-2), so it is
        // charged as its own squarings and multiplications.
        let mul_cost = c.mul * guest_cost::MUL;
        let sq_cost = c.sq * guest_cost::SQ;
        let add_cost = c.add * guest_cost::ADD;
        let inv_cost = c.inv * guest_cost::INV;
        let arith_cost = mul_cost + sq_cost + add_cost + inv_cost;
        let total = hash_cost + arith_cost;
        let pct = |x: u64| 100.0 * x as f64 / total as f64;

        println!("== {NAME}, {num_signatures} signature(s), num_bits={num_bits}, verify_ok={} ==", res.is_ok());
        println!("  {:<22} {:>12} {:>7} {:>14} {:>7}", "operation", "count", "instr", "instructions", "share");
        for (label, count, unit, cost) in [
            ("poseidon2 perm", perms, guest_cost::PERM, perm_cost),
            ("limb split (hash glue)", splits, guest_cost::SPLIT, split_cost),
            ("ed25519 mul", c.mul, guest_cost::MUL, mul_cost),
            ("ed25519 square", c.sq, guest_cost::SQ, sq_cost),
            ("ed25519 add/sub", c.add, guest_cost::ADD, add_cost),
            ("ed25519 inverse", c.inv, guest_cost::INV, inv_cost),
        ] {
            println!("  {label:<22} {count:>12} {unit:>7} {cost:>14} {:>6.1}%", pct(cost));
        }
        println!("  {:-<66}", "");
        println!("  {:<22} {:>12} {:>7} {:>14} {:>6.1}%", "HASHING", "", "", hash_cost, pct(hash_cost));
        println!("  {:<22} {:>12} {:>7} {:>14} {:>6.1}%", "ARITHMETIC", "", "", arith_cost, pct(arith_cost));
        println!("  {:<22} {:>12} {:>7} {:>14}", "modelled total", "", "", total);
    }
    #[cfg(not(feature = "count-ops"))]
    println!("{NAME} verify_ok={:?} (build with --features count-ops)", res.is_ok());
}
