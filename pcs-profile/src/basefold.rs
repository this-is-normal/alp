//! Basefold PCS-only profile: drives the core `BasefoldPcs::{commit, open, verify_opening}` on a
//! random multilinear at a random point (rate 1/4 foldable code, BabyBear Merkle hash).

use p3_basefold::{BasefoldPcs, RandomFoldableCode, eval_multilinear};
use p3_ed25519::Ed25519;
use p3_stark_stats::{HashCounts, serialized_size};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::{PcsMetrics, ProfileablePcs, build_challenger, build_mmcs, reset_arith, snapshot_arith};

/// Basefold over a rate-`1/inv_rate` random foldable code.
pub struct BasefoldProfile {
    pub inv_rate: usize,
    pub code_seed: u64,
}

impl Default for BasefoldProfile {
    fn default() -> Self {
        Self { inv_rate: 4, code_seed: 0xBA5E_F01D }
    }
}

impl ProfileablePcs for BasefoldProfile {
    fn name(&self) -> &'static str {
        "Basefold"
    }

    fn run(&self, num_vars: usize, q: usize, seed: u64) -> PcsMetrics {
        let counts = HashCounts::new();
        let mmcs = build_mmcs(&counts);
        let code = RandomFoldableCode::<Ed25519>::new(num_vars, self.inv_rate, self.code_seed);
        let qc = q.min(self.inv_rate << num_vars);
        // Rate-1/4 foldable code, list-decoding proximity radius ≈ 0.215.
        let sec_bits = qc as f64 * -(1.0 - 0.215_f64).log2();
        let pcs = BasefoldPcs::new(code, mmcs, qc);

        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let evals: Vec<Ed25519> = (0..(1usize << num_vars)).map(|_| rng.random()).collect();
        let z: Vec<Ed25519> = (0..num_vars).map(|_| rng.random()).collect();
        // NB: `open` returns the folded final value, not f(z); the standalone verify wants f(z).
        let y = eval_multilinear(&evals, &z);

        let mut ch_p = build_challenger(&counts);
        let t_commit = std::time::Instant::now();
        let (comm, pd) = pcs.commit(evals, &mut ch_p);
        let commit_us = t_commit.elapsed().as_micros();
        let t_open = std::time::Instant::now();
        let (_folded, proof) = pcs.open(&pd, &z, &mut ch_p);
        let open_us = t_open.elapsed().as_micros();
        let proof_bytes = serialized_size(&proof);

        // Measure the VERIFY only.
        counts.reset();
        reset_arith();
        let mut ch_v = build_challenger(&counts);
        let t_verify = std::time::Instant::now();
        pcs.verify(&comm, &z, y, &proof, &mut ch_v).expect("basefold verify");
        let verify_us = t_verify.elapsed().as_micros();

        let (hl, hc, hs) = counts.snapshot();
        let (mul, sq, inv, add) = snapshot_arith();
        PcsMetrics {
            pcs: "Basefold",
            regime: "list",
            num_vars,
            q: qc,
            sec_bits,
            proof_bytes,
            hash_leaf: hl as u64,
            hash_compress: hc as u64,
            hash_sponge: hs as u64,
            mul,
            sq,
            inv,
            add,
            commit_us,
            open_us,
            verify_us,
        }
    }
}
