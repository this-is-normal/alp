//! Brakedown tensor-PCS-only profile: drives the core `TensorPcs::{commit, open, verify}` on a
//! random multilinear at a random point, using a vetted [`Preset`] (small-proof / balanced /
//! fast-prover) and a real **distance-derived** query count — `t = ⌈λ/−log₂(1−e)⌉`, e = δ/3
//! (provable) or δ/2 (conjectured) — clamped to the codeword length.

use p3_brakedown::{BrakedownCode, LinearCode, Preset, TensorPcs};
use p3_ed25519::Ed25519;
use p3_multilinear_util::point::Point;
use p3_stark_stats::{HashCounts, serialized_size};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::{Mmcs, PcsMetrics, ProfileablePcs, build_challenger, build_mmcs, reset_arith, snapshot_arith};

/// Brakedown tensor PCS at one operating point + query regime.
pub struct BrakedownProfile {
    pub preset: Preset,
    /// `false` = provable δ/3; `true` = conjectured δ/2 (Diamond–Posen Conjecture-1).
    pub conjectured: bool,
    pub security_bits: usize,
    pub code_seed: u64,
}

impl BrakedownProfile {
    #[must_use]
    pub fn new(preset: Preset, conjectured: bool) -> Self {
        Self { preset, conjectured, security_bits: 100, code_seed: 0xB2A4E_D0 }
    }
}

impl ProfileablePcs for BrakedownProfile {
    fn name(&self) -> &'static str {
        self.preset.name
    }

    fn run(&self, num_vars: usize, _q_ignored: usize, seed: u64) -> PcsMetrics {
        let counts = HashCounts::new();
        let mmcs = build_mmcs(&counts);

        // Real distance-derived query count, then the paper's optimal row/column split for it.
        let t = self.preset.num_queries(self.security_bits, self.conjectured);
        let log_rows = p3_brakedown::params::optimal_log_rows(num_vars, t, &self.preset.params);
        let cols = 1usize << (num_vars - log_rows);
        let mut code_rng = ChaCha8Rng::seed_from_u64(self.code_seed);
        let code = BrakedownCode::<Ed25519, Ed25519>::from_message_len_with_params(
            cols,
            &self.preset.params,
            &mut code_rng,
        );
        let qc = t.min(code.codeword_len());
        // Actual soundness with the (possibly clamped) query count.
        let e = if self.conjectured { self.preset.delta / 2.0 } else { self.preset.delta / 3.0 };
        let sec_bits = qc as f64 * -(1.0 - e).log2();
        let pcs = TensorPcs::<Ed25519, Ed25519, BrakedownCode<Ed25519, Ed25519>, Mmcs>::new(
            code, mmcs, log_rows, qc,
        );

        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let evals: Vec<Ed25519> = (0..(1usize << num_vars)).map(|_| rng.random()).collect();
        let point = Point::new((0..num_vars).map(|_| rng.random::<Ed25519>()).collect());

        let mut ch_p = build_challenger(&counts);
        let t_commit = std::time::Instant::now();
        let (comm, pd) = pcs.commit(evals, &mut ch_p);
        let commit_us = t_commit.elapsed().as_micros();
        let t_open = std::time::Instant::now();
        let (value, proof) = pcs.open(&pd, &point, &mut ch_p);
        let open_us = t_open.elapsed().as_micros();
        let proof_bytes = serialized_size(&proof);

        // Measure the VERIFY only.
        counts.reset();
        reset_arith();
        let mut ch_v = build_challenger(&counts);
        let t_verify = std::time::Instant::now();
        pcs.verify(&comm, &point, value, &proof, &mut ch_v).expect("brakedown verify");
        let verify_us = t_verify.elapsed().as_micros();

        let (hl, hc, hs) = counts.snapshot();
        let (mul, sq, inv, add) = snapshot_arith();
        PcsMetrics {
            pcs: self.preset.name,
            regime: if self.conjectured { "d/2" } else { "d/3" },
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
