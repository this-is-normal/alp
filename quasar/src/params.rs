//! Shape and security parameters.

use p3_security::QaRegime;
use serde::{Deserialize, Serialize};

/// The public shape of a Quasar instance.
///
/// The committed polynomial has `2^{log_rows + log_cols}` evaluations, arranged
/// as `N₁ = 2^log_rows` rows of `N₂ = 2^log_cols`. Each row is QA-encoded to
/// `index · N₂` symbols and the columns are Merkle-committed.
///
/// `log_rows` is the central trade-off: the proof carries `num_queries · N₁`
/// field elements from the opened columns, so small `N₁` means small proofs —
/// but the code-switching layer then works on longer rows (`N₂ = N/N₁`), which
/// is where the linear prover cost sits. The paper uses `N₁ = 2⁸`; its own
/// benchmarks use `2⁶`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct QuasarParams {
    /// `log₂N₁` — matrix rows, i.e. the height of an opened column.
    pub log_rows: usize,
    /// `log₂N₂` — matrix columns, i.e. the QA message length.
    pub log_cols: usize,
    /// QA index `c` = inverse rate.
    pub index: usize,
    /// Committed-matrix columns opened in the proximity/consistency test.
    pub num_queries: usize,
}

impl QuasarParams {
    /// Build from an explicit shape.
    pub const fn new(log_rows: usize, log_cols: usize, index: usize, num_queries: usize) -> Self {
        Self {
            log_rows,
            log_cols,
            index,
            num_queries,
        }
    }

    /// Build with `num_queries` derived from the proven QA distance.
    ///
    /// # Returns
    ///
    /// `None` when the field is too small for a useful QA distance at this row
    /// length (any 31-bit field, for instance).
    pub fn from_security(
        log_rows: usize,
        log_cols: usize,
        index: usize,
        field_bits: usize,
        security_bits: usize,
        distance_failure_bits: usize,
    ) -> Option<Self> {
        let regime = QaRegime {
            row_log_size: log_cols,
            index,
            field_bits,
            distance_failure_bits,
            query_security_bits: security_bits,
        };
        Some(Self::new(log_rows, log_cols, index, regime.num_queries()?))
    }

    /// Choose the row/column split that minimises proof size for a committed
    /// polynomial of `num_vars` variables.
    ///
    /// The proof has two size terms that pull in opposite directions:
    ///
    /// ```text
    ///   t · N₁                 the opened columns          (grows with N₁)
    ///   O(q_bf · log²N₂)       the two BaseFold openings   (grows with N₂)
    /// ```
    ///
    /// so a *narrow* matrix keeps the columns short while a *short* one keeps
    /// the code-switching openings cheap. This scores both in bytes, over every
    /// power-of-two aspect ratio whose codeword can actually hold `t` columns.
    ///
    /// The score mirrors the `field_elems`/`digests` terms of
    /// `quasar_cost_model.py`, whose predictions are within 1 % of measured
    /// proof bytes; `chooser_matches_cost_model` pins the two together.
    ///
    /// # Returns
    ///
    /// `None` when the field admits no useful QA distance, or when no split
    /// leaves at least one column variable.
    pub fn choose(
        num_vars: usize,
        index: usize,
        field_bits: usize,
        security_bits: usize,
        distance_failure_bits: usize,
        basefold_inv_rate: usize,
        basefold_queries: usize,
    ) -> Option<Self> {
        Self::search(
            num_vars,
            index,
            field_bits,
            security_bits,
            distance_failure_bits,
            basefold_inv_rate,
            basefold_queries,
            true,
        )
    }

    /// [`Self::choose`] without the security filter: the smallest-proof shape at
    /// whatever soundness the codeword can carry, with `num_queries` clamped.
    ///
    /// For demo-sized instances, where no shape reaches the target and the
    /// alternative is not proving at all. Callers must report
    /// [`Self::achieved_security_bits`] rather than implying the target.
    ///
    /// # Returns
    ///
    /// `None` only if the field admits no useful QA distance at all.
    pub fn choose_best_effort(
        num_vars: usize,
        index: usize,
        field_bits: usize,
        security_bits: usize,
        distance_failure_bits: usize,
        basefold_inv_rate: usize,
        basefold_queries: usize,
    ) -> Option<Self> {
        let params = Self::search(
            num_vars,
            index,
            field_bits,
            security_bits,
            distance_failure_bits,
            basefold_inv_rate,
            basefold_queries,
            false,
        )?;
        Some(Self::new(
            params.log_rows,
            params.log_cols,
            index,
            params.num_queries.min(params.codeword_len()),
        ))
    }

    #[allow(clippy::too_many_arguments)]
    fn search(
        num_vars: usize,
        index: usize,
        field_bits: usize,
        security_bits: usize,
        distance_failure_bits: usize,
        basefold_inv_rate: usize,
        basefold_queries: usize,
        require_target: bool,
    ) -> Option<Self> {
        let mut best: Option<(u128, Self)> = None;
        for log_rows in 0..num_vars {
            let log_cols = num_vars - log_rows;
            let candidate = Self::from_security(
                log_rows,
                log_cols,
                index,
                field_bits,
                security_bits,
                distance_failure_bits,
            )?;
            // A codeword shorter than the query count cannot reach the target.
            if require_target && candidate.num_queries > candidate.codeword_len() {
                continue;
            }
            let clamped = Self::new(
                log_rows,
                log_cols,
                index,
                candidate.num_queries.min(candidate.codeword_len()),
            );
            let cost = clamped.proof_size_score(basefold_inv_rate, basefold_queries);
            if best.as_ref().is_none_or(|(b, _)| cost < *b) {
                best = Some((cost, candidate));
            }
        }
        best.map(|(_, params)| params)
    }

    /// Predicted serialized proof size in bytes, used by [`Self::choose`].
    ///
    /// Weights match the validated cost model: 33 bytes per field element
    /// (postcard writes a length byte then 32) and 39 per digest (8 varint
    /// BabyBear limbs, measured mean). A pure-ed25519 hash has 33-byte digests
    /// instead; across every size tested that moves the argmin in 2 cases out
    /// of 16, and by one step worth under 0.5% of proof size — the optimum is
    /// flat, so one weighting serves both hash stacks.
    pub fn proof_size_score(&self, basefold_inv_rate: usize, basefold_queries: usize) -> u128 {
        const FIELD_BYTES: u128 = 33;
        const DIGEST_BYTES: u128 = 39;

        let (n2, c, t) = (self.log_cols, self.index, self.num_queries as u128);
        let d = n2 as u128; // BaseFold fold rounds
        let log_l = self.log_codeword_len() as u128;
        let log_n0 = (basefold_inv_rate << n2).trailing_zeros() as u128;
        let qbf = basefold_queries.min(basefold_inv_rate << n2) as u128;
        // Per BaseFold query: paths of shrinking depth over the folded oracles.
        let path_sum = d * log_n0 - d * (d + 1) / 2;

        // Field elements: the opened columns, each opening's fold sumcheck and
        // opened rows, and the O(log N) sumcheck transcripts and claims.
        let opened_rows = |k: u128| 3 * d + 1 + qbf * (2 * k + 2 * (d - 1));
        let elements = t * self.num_rows() as u128
            + opened_rows(2 * c as u128)
            + opened_rows(c as u128 - 1)
            + 4
            + 2 * log_l
            + 7 * d
            + 10 * c as u128;
        // Digests: one path per opened column, plus each opening's fold roots
        // and query paths.
        let digests = t * log_l + 2 * ((d - 1) + qbf * path_sum);

        elements * FIELD_BYTES + digests * DIGEST_BYTES
    }

    /// Number of variables of the committed multilinear polynomial.
    pub const fn num_vars(&self) -> usize {
        self.log_rows + self.log_cols
    }

    /// `N₁`, the number of matrix rows.
    pub const fn num_rows(&self) -> usize {
        1 << self.log_rows
    }

    /// `N₂`, the number of matrix columns (the QA message length).
    pub const fn num_cols(&self) -> usize {
        1 << self.log_cols
    }

    /// `c·N₂`, the length of one encoded row — and the number of committed
    /// columns.
    pub const fn codeword_len(&self) -> usize {
        self.index << self.log_cols
    }

    /// `log₂(c·N₂)`, the variable count of the concatenated-codeword
    /// polynomial the selector sumcheck runs over.
    ///
    /// # Panics
    ///
    /// Panics unless `index` is a power of two, which the block-indexing
    /// variables require.
    pub const fn log_codeword_len(&self) -> usize {
        assert!(
            self.index.is_power_of_two(),
            "QA index must be a power of two"
        );
        self.log_cols + self.index.trailing_zeros() as usize
    }

    /// `log₂c` — the number of leading variables that select a codeword block.
    pub const fn log_index(&self) -> usize {
        self.index.trailing_zeros() as usize
    }

    /// Number of length-`N₂` blocks committed online: the `c` blocks of the
    /// proximity codeword plus the `c` of the evaluation codeword.
    pub const fn num_online_blocks(&self) -> usize {
        2 * self.index
    }

    /// Number of preprocessed coefficient vectors, `c − 1`.
    pub const fn num_coeff_vectors(&self) -> usize {
        self.index - 1
    }

    /// Bits of column-sampling soundness this shape **actually** delivers:
    /// `−t·log₂(1 − δ/3)`.
    ///
    /// This is not the target. Small instances clamp `num_queries` to the
    /// codeword length, and the achieved figure then falls far below it — a
    /// 2⁷-element codeword can carry at most 128 of the ~420 columns a 100-bit
    /// instance needs. Configurations should report this rather than imply the
    /// target, and tests at demo sizes should expect it to be low.
    ///
    /// # Returns
    ///
    /// `None` if the field admits no useful QA distance.
    pub fn achieved_security_bits(&self, field_bits: usize, distance_failure_bits: usize) -> Option<f64> {
        let regime = QaRegime {
            row_log_size: self.log_cols,
            index: self.index,
            field_bits,
            distance_failure_bits,
            query_security_bits: 0,
        };
        regime
            .sampling_error(self.num_queries)
            .map(|bits| bits.bits())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shape_arithmetic() {
        let p = QuasarParams::new(6, 14, 2, 425);
        assert_eq!(p.num_vars(), 20);
        assert_eq!(p.num_rows(), 64);
        assert_eq!(p.num_cols(), 1 << 14);
        assert_eq!(p.codeword_len(), 1 << 15);
        assert_eq!(p.log_codeword_len(), 15);
        assert_eq!(p.log_index(), 1);
        assert_eq!(p.num_online_blocks(), 4);
        assert_eq!(p.num_coeff_vectors(), 1);
    }

    /// The chooser must land on the same shape a brute-force sweep of
    /// `quasar_cost_model.py` picks by predicted proof *bytes* — the score here
    /// is only a proxy (digests counted as field elements), so agreeing on the
    /// argmin is the property that matters.
    ///
    /// Reference values, ed25519 at the paper's targets, BaseFold rate 1/4 with
    /// 288 repetitions:
    /// `(num_vars, index) -> log_rows`.
    #[test]
    fn chooser_matches_cost_model() {
        const REFERENCE: [(usize, usize, usize); 12] = [
            (14, 2, 5),
            (14, 4, 6),
            (16, 2, 5),
            (16, 4, 6),
            (18, 2, 5),
            (18, 4, 6),
            (20, 2, 5),
            (20, 4, 6),
            (25, 2, 6),
            (25, 4, 6),
            (29, 2, 6),
            (29, 4, 7),
        ];
        for (num_vars, index, expected) in REFERENCE {
            let chosen = QuasarParams::choose(num_vars, index, 255, 103, 103, 4, 288)
                .unwrap_or_else(|| panic!("no shape for 2^{num_vars} at index {index}"));
            assert_eq!(
                chosen.log_rows, expected,
                "2^{num_vars} index {index}: chose N1=2^{}, cost model says 2^{expected}",
                chosen.log_rows
            );
            assert_eq!(chosen.num_vars(), num_vars);
            // The whole point of the filter: the codeword must be able to hold
            // the query count the target demands.
            assert!(chosen.num_queries <= chosen.codeword_len());
        }
    }

    /// Below 2⁸ evaluations no shape can hold the ~410 columns a 100-bit
    /// instance needs, so the chooser declines rather than quietly returning an
    /// under-secured shape.
    #[test]
    fn chooser_declines_when_no_shape_reaches_the_target() {
        // A rate-1/2 codeword needs N₂ >= t/2 ~ 205 columns, i.e. n₂ >= 8.
        for num_vars in 1..8 {
            assert_eq!(
                QuasarParams::choose(num_vars, 2, 255, 103, 103, 4, 288),
                None,
                "2^{num_vars} should admit no fully secure shape"
            );
        }
        let smallest = QuasarParams::choose(8, 2, 255, 103, 103, 4, 288)
            .expect("2^8 is the smallest fully secure size at index 2");
        assert_eq!((smallest.log_rows, smallest.log_cols), (0, 8));

        // Best-effort still answers below the threshold, with a clamped count.
        let demo = QuasarParams::choose_best_effort(7, 2, 255, 103, 103, 4, 288).unwrap();
        assert_eq!(demo.num_vars(), 7);
        assert_eq!(demo.num_queries, demo.num_queries.min(demo.codeword_len()));
        assert!(demo.achieved_security_bits(255, 103).unwrap() < 103.0);
    }

    /// A clamped query count must report the soundness it actually delivers,
    /// not the target it was asked for.
    #[test]
    fn achieved_security_is_honest_about_clamping() {
        // A real instance reaches its target.
        let real = QuasarParams::from_security(6, 16, 2, 255, 103, 103).unwrap();
        let bits = real.achieved_security_bits(255, 103).unwrap();
        assert!((bits - 103.0).abs() < 1.0, "got {bits}");

        // A demo-sized one, clamped to a 2⁷-element codeword, does not.
        let demo = QuasarParams::new(3, 6, 2, 128);
        let bits = demo.achieved_security_bits(255, 103).unwrap();
        assert!(bits < 40.0, "clamped instance should be far below target, got {bits}");
    }

    #[test]
    fn security_derived_queries_match_the_regime() {
        // ed25519 at the paper's targets.
        let p = QuasarParams::from_security(6, 16, 2, 255, 103, 103).unwrap();
        assert_eq!(p.num_queries, 425);
        // Rate 1/4 needs far fewer columns.
        let p4 = QuasarParams::from_security(6, 16, 4, 255, 103, 103).unwrap();
        assert_eq!(p4.num_queries, 256);
        // A 31-bit field has no admissible distance, so no parameters.
        assert!(QuasarParams::from_security(6, 16, 2, 31, 100, 100).is_none());
    }
}
