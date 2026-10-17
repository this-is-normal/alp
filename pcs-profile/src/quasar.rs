//! Quasar PCS-only profile: drives `QuasarPcs::{commit, open, verify}` on a
//! random multilinear at a random point, under the same BabyBear-hash stack the
//! Basefold and Brakedown profiles use.
//!
//! Two query counts matter here and are reported separately from the `q`
//! column: `q` is the number of committed *columns* opened (derived from the QA
//! distance bound, `p3-security::QaRegime`), while the code-switching openings
//! use `bf_queries` BaseFold repetitions.

use p3_basefold::{BasefoldPcs, RandomFoldableCode};
use p3_ed25519::Ed25519;
use p3_qa_code::QaCode;
use p3_quasar::{QuasarParams, QuasarPcs, eval_mle};
use p3_security::QaRegime;
use p3_stark_stats::{HashCounts, serialized_size};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::{PcsMetrics, ProfileablePcs, build_challenger, build_mmcs, reset_arith, snapshot_arith};

/// Quasar over a QA code of the given inverse rate, code switching into a
/// rate-`1/bf_inv_rate` foldable code.
pub struct QuasarProfile {
    /// `log₂N₁` — matrix rows, i.e. the height of an opened column.
    pub log_rows: usize,
    /// QA index `c` (inverse rate) of the committing code.
    pub index: usize,
    /// BaseFold inverse rate for the code-switching openings.
    pub bf_inv_rate: usize,
    /// BaseFold query repetitions.
    pub bf_queries: usize,
    pub code_seed: u64,
}

impl Default for QuasarProfile {
    fn default() -> Self {
        Self {
            log_rows: 6,
            index: 2,
            bf_inv_rate: 4,
            bf_queries: 288,
            code_seed: 0x0A_C0DE,
        }
    }
}

impl QuasarProfile {
    /// The column count the proven QA distance demands at this row length over
    /// ed25519, at the paper's targets (λ_samp = λ_dist = 103).
    ///
    /// # Returns
    ///
    /// `None` if the parameters admit no useful QA distance — which cannot
    /// happen over a 255-bit field, but is checked rather than assumed.
    pub fn security_queries(&self, num_vars: usize) -> Option<usize> {
        QaRegime {
            row_log_size: num_vars.saturating_sub(self.log_rows),
            index: self.index,
            field_bits: 255,
            distance_failure_bits: 103,
            query_security_bits: 103,
        }
        .num_queries()
    }
}

impl ProfileablePcs for QuasarProfile {
    fn name(&self) -> &'static str {
        "Quasar"
    }

    fn run(&self, num_vars: usize, q: usize, seed: u64) -> PcsMetrics {
        assert!(
            num_vars > self.log_rows,
            "need at least two columns: num_vars must exceed log_rows"
        );
        let log_cols = num_vars - self.log_rows;

        let counts = HashCounts::new();
        // The column tree and the Basefold trees share the hash stack, so the
        // perm counters cover both.
        let qc = q.min(self.index << log_cols);
        let params = QuasarParams::new(self.log_rows, log_cols, self.index, qc);
        let code = QaCode::<Ed25519>::new_random(
            log_cols,
            self.index,
            &mut ChaCha8Rng::seed_from_u64(self.code_seed),
        );
        let basefold = BasefoldPcs::new(
            RandomFoldableCode::<Ed25519>::new(log_cols, self.bf_inv_rate, 0xBA5E_F01D),
            build_mmcs(&counts),
            self.bf_queries.min(self.bf_inv_rate << log_cols),
        );
        let (pcs, pk) = QuasarPcs::new(params, code, build_mmcs(&counts), basefold);

        // Provable column-sampling regime: each opened column contributes
        // −log₂(1 − δ/3) bits.
        let delta = QaRegime {
            row_log_size: log_cols,
            index: self.index,
            field_bits: 255,
            distance_failure_bits: 103,
            query_security_bits: 103,
        }
        .distance()
        .expect("ed25519 admits a QA distance bound");
        let sec_bits = qc as f64 * -(1.0 - delta / 3.0).log2();

        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let evals: Vec<Ed25519> = (0..(1usize << num_vars)).map(|_| rng.random()).collect();
        let z: Vec<Ed25519> = (0..num_vars).map(|_| rng.random()).collect();
        let expected = eval_mle(&evals, &z);

        let mut ch_p = build_challenger(&counts);
        let t_commit = std::time::Instant::now();
        let (comm, pd) = pcs.commit(evals, &mut ch_p);
        let commit_us = t_commit.elapsed().as_micros();

        let t_open = std::time::Instant::now();
        let (value, proof) = pcs.open(&pk, &pd, &z, &mut ch_p);
        let open_us = t_open.elapsed().as_micros();
        assert_eq!(value, expected, "Quasar opened the wrong value");
        let proof_bytes = serialized_size(&proof);

        // Measure the VERIFY only.
        counts.reset();
        reset_arith();
        let mut ch_v = build_challenger(&counts);
        let t_verify = std::time::Instant::now();
        pcs.verify(&comm, &z, value, &proof, &mut ch_v)
            .expect("quasar verify");
        let verify_us = t_verify.elapsed().as_micros();

        let (hl, hc, hs) = counts.snapshot();
        let (mul, sq, inv, add) = snapshot_arith();
        PcsMetrics {
            pcs: "Quasar",
            regime: "d/3",
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
