//! Quasi-abelian (QA) code soundness: the distance bound for a *random*
//! rank-1 systematic QA code, and the column-query count it feeds.
//!
//! A rank-1 systematic QA code of index `c` and message length `N = 2ⁿ`
//! has generator `(I | A_2 | … | A_c)`, where each `A_i` is drawn from the
//! recursive dyadic distribution
//!
//! ```text
//!   A_0 = [ξ],   ξ ← F                    A_ℓ = | A_{ℓ-1,L}  A_{ℓ-1,R} |
//!                                               | A_{ℓ-1,R}  A_{ℓ-1,L} |
//! ```
//!
//! with the two halves sampled independently. Its relative distance is a
//! random variable; the bound below is the tail estimate of [LLX⁺26b]
//! (Thm 4.4 / Cor 3.21), as used by the Quasar paper. Note that it is
//! *field-size hungry*: the `(1 + log n / log p) / c` term makes the bound
//! vacuous for 31-bit fields at any useful length, marginal at 64 bits, and
//! comfortable from 127 bits up.
//!
//! The committed-matrix column test opens `t` columns and needs
//! `(1 − δ/3)ᵗ ≤ 2^-λ`; the `δ/3` testing radius is the provable regime
//! (there is no list-decoding conjecture in play here).
//!
//! # References
//! - Quasar: field-agnostic PCS from quasi-abelian codes
//!   ([2026/1839](https://eprint.iacr.org/2026/1839)) — §2, App. C
//! - Li–Li–Xing⁺ QA codes, referenced there as [LLX⁺26b] — Thm 4.4, Cor 3.21

use libm::{ceil, expm1, log, log2, pow};
use serde::Serialize;

use crate::error::ErrorBits;

/// Lower edge of the bisection search for `δ`; a returned distance at or
/// below this means "no admissible distance for these parameters".
const DELTA_FLOOR: f64 = 1e-6;

/// Bisection steps used by [`distance_lower_bound`]. 100 halvings of
/// `(0, 1)` is far below `f64` resolution, so the result is exact to
/// double precision.
const BISECTION_STEPS: usize = 100;

/// The `g_p(δ)` term of the QA distance bound:
/// `1 − δ + (δ·log₂δ + (1−δ)·log₂(1−δ)) / log₂|F|`.
///
/// This is the entropy-corrected agreement exponent; the bound is only
/// meaningful while it exceeds `(1 + log₂N / log₂|F|) / c`.
#[must_use]
pub fn gp(delta: f64, field_bits: usize) -> f64 {
    debug_assert!(delta > 0.0 && delta < 1.0);
    let bits = field_bits as f64;
    1.0 - delta + (delta * log2(delta) + (1.0 - delta) * log2(1.0 - delta)) / bits
}

/// `log₂ Pr[δ(C_QA) < δ]` for a random rank-1 systematic QA code of index
/// `c` and message length `2^row_log_size` over a `field_bits`-bit field.
///
/// Returns `+∞` when the bound is vacuous (the `ε > 0` precondition fails),
/// which is the honest answer: the estimate says nothing at those
/// parameters. Takes the tighter of the two union-bound branches.
#[must_use]
pub fn distance_failure_log2(
    delta: f64,
    row_log_size: usize,
    index: usize,
    field_bits: usize,
) -> f64 {
    assert!(index >= 2, "QA index (inverse rate) must be at least 2");
    if !(delta > 0.0 && delta < 1.0) {
        return f64::INFINITY;
    }

    let c = index as f64;
    let log_n = row_log_size as f64;
    let log_p = field_bits as f64;

    let eps = gp(delta, field_bits) - (1.0 + log_n / log_p) / c;
    if eps <= 0.0 {
        return f64::INFINITY;
    }

    // log2(1 - 2^{-log_p·c·ε}), the geometric-series denominator. Saturates
    // to 0 once the exponent is large enough that the term is 1.
    let x = log_p * c * eps;
    let denom_log = if x < 60.0 {
        log2(-expm1(-log(2.0) * x))
    } else {
        0.0
    };

    // Branch 1: pairs of parity blocks.
    let t1a = log2(c * (c - 1.0) / 2.0) + log_n - 2.0 * log_p;
    let t1b = -log_p * ceil((c - 1.0) / (c * delta)) * c * eps - denom_log - log_p;

    // Branch 2: single parity blocks.
    let t2a = log2(c) + log_n - log_p;
    let t2b = -log_p * ceil(1.0 / delta) * c * eps - denom_log - log_p;

    log2_add(t1a, t1b).min(log2_add(t2a, t2b))
}

/// `log₂(2^a + 2^b)`, computed without intermediate overflow.
fn log2_add(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    if !m.is_finite() {
        return m;
    }
    m + log2(pow(2.0, a - m) + pow(2.0, b - m))
}

/// The largest relative distance `δ` whose failure probability clears
/// `failure_bits` bits, i.e. `Pr[δ(C_QA) < δ] ≤ 2^-failure_bits`.
///
/// # Returns
///
/// - `Some(δ)` for parameters admitting a useful bound.
/// - `None` when the field is too small for the length (e.g. any 31-bit
///   field), rather than a meaninglessly tiny `δ`.
#[must_use]
pub fn distance_lower_bound(
    row_log_size: usize,
    index: usize,
    field_bits: usize,
    failure_bits: usize,
) -> Option<f64> {
    let target = -(failure_bits as f64);
    let mut lo = 1e-9_f64;
    let mut hi = 1.0 - 1e-9_f64;

    for _ in 0..BISECTION_STEPS {
        let mid = 0.5 * (lo + hi);
        if distance_failure_log2(mid, row_log_size, index, field_bits) <= target {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    (lo > DELTA_FLOOR).then_some(lo)
}

/// Columns to open so that `(1 − δ/3)ᵗ ≤ 2^-security_bits`:
/// `t = ⌈security_bits / −log₂(1 − δ/3)⌉`.
#[must_use]
pub fn num_queries_from_distance(delta: f64, security_bits: usize) -> usize {
    assert!(delta > 0.0 && delta < 1.0, "distance must lie in (0, 1)");
    let per_query = -log2(1.0 - delta / 3.0);
    ceil(security_bits as f64 / per_query) as usize
}

/// The security-relevant shape of a QA-code commitment.
///
/// `row_log_size` is `log₂N₂` — the *row* length of the committed matrix,
/// which is what the code is applied to; the number of matrix rows `N₁`
/// does not enter the distance bound (it only scales proof size).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize)]
pub struct QaRegime {
    /// `log₂` of the QA message length (one matrix row).
    pub row_log_size: usize,
    /// QA index `c` = inverse rate; the codeword is `c · 2^row_log_size`.
    pub index: usize,
    /// `log₂|F|` of the field the code is over.
    pub field_bits: usize,
    /// Bits of failure probability allowed for the random-code distance
    /// bound itself (`ε_dist`).
    pub distance_failure_bits: usize,
    /// Target bits for the column-sampling error (`ε_samp`), which fixes
    /// the query count.
    pub query_security_bits: usize,
}

impl QaRegime {
    /// Codeword length `c · 2^row_log_size`, in field elements.
    #[must_use]
    pub const fn codeword_len(&self) -> usize {
        self.index << self.row_log_size
    }

    /// The proven relative distance, or `None` if the parameters admit no
    /// useful bound. See [`distance_lower_bound`].
    #[must_use]
    pub fn distance(&self) -> Option<f64> {
        distance_lower_bound(
            self.row_log_size,
            self.index,
            self.field_bits,
            self.distance_failure_bits,
        )
    }

    /// Number of committed-matrix columns to open.
    #[must_use]
    pub fn num_queries(&self) -> Option<usize> {
        self.distance()
            .map(|d| num_queries_from_distance(d, self.query_security_bits))
    }

    /// `ε_dist` — the distance bound's own failure probability, evaluated
    /// at the `δ` this regime uses.
    #[must_use]
    pub fn distance_error(&self) -> Option<ErrorBits> {
        let d = self.distance()?;
        Some(ErrorBits::from_log2(-distance_failure_log2(
            d,
            self.row_log_size,
            self.index,
            self.field_bits,
        )))
    }

    /// `ε_samp = (1 − δ/3)ᵗ` — a far-from-codeword matrix surviving every
    /// sampled column.
    #[must_use]
    pub fn sampling_error(&self, num_queries: usize) -> Option<ErrorBits> {
        let d = self.distance()?;
        Some(ErrorBits::from_log2(
            -(num_queries as f64) * log2(1.0 - d / 3.0),
        ))
    }

    /// `ε_aux = 2·(1 − 2δ/3)ᵗ` — the auxiliary sampled-column term of the
    /// Quasar analysis (App. C, Eq. 2).
    #[must_use]
    pub fn auxiliary_error(&self, num_queries: usize) -> Option<ErrorBits> {
        let d = self.distance()?;
        Some(ErrorBits::from_log2(
            -1.0 - (num_queries as f64) * log2(1.0 - 2.0 * d / 3.0),
        ))
    }

    /// `ε_fld` — the collected `1/|F|` terms (sumcheck rounds, batching
    /// challenges, selector weights, extractor linear algebra), per Quasar
    /// App. C:
    /// `(12n₂ + 2c − 4 + (t − 1 + 4ℓ) + (δL/3 + 1)) / |F|`, `ℓ = ⌈log₂L⌉`.
    #[must_use]
    pub fn field_error(&self, num_queries: usize) -> Option<ErrorBits> {
        let d = self.distance()?;
        let l = self.codeword_len() as f64;
        let log_l = ceil(log2(l));
        let numerator = 12.0 * self.row_log_size as f64
            + 2.0 * self.index as f64
            - 4.0
            + (num_queries as f64 - 1.0 + 4.0 * log_l)
            + (d * l / 3.0 + 1.0);
        Some(ErrorBits::from_log2(
            self.field_bits as f64 - log2(numerator),
        ))
    }

    /// The complete Quasar union bound (App. C, Eq. 2):
    /// `ε_dist + ε_samp + ε_aux + ε_fld + ε_BF + ε_VC`.
    ///
    /// `basefold_error` and `vector_commitment_error` come from the
    /// code-switching target and the Merkle commitment respectively — they
    /// are not QA-code properties, so the caller supplies them.
    #[must_use]
    pub fn total_error(
        &self,
        basefold_error: ErrorBits,
        vector_commitment_error: ErrorBits,
    ) -> Option<ErrorBits> {
        let t = self.num_queries()?;
        Some(ErrorBits::sum(&[
            self.distance_error()?,
            self.sampling_error(t)?,
            self.auxiliary_error(t)?,
            self.field_error(t)?,
            basefold_error,
            vector_commitment_error,
        ]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Paper Table 4: `(k, n₂, δ_QA, t, λ_fld, λ_total)` over Mersenne127
    /// with `c = 2`, `λ_dist = λ_samp target = 103`, 249 BaseFold
    /// repetitions and a 128-bit vector commitment.
    const TABLE_4: [(usize, usize, f64, usize, f64, f64); 10] = [
        (20, 12, 0.444951, 445, 116.1, 101.5),
        (21, 13, 0.441024, 450, 115.4, 101.6),
        (22, 14, 0.437098, 454, 114.6, 101.6),
        (23, 15, 0.433173, 458, 113.7, 101.5),
        (24, 16, 0.429248, 463, 112.8, 101.6),
        (25, 17, 0.425324, 467, 111.8, 101.5),
        (26, 18, 0.421401, 472, 110.8, 101.5),
        (27, 19, 0.417478, 477, 109.8, 101.6),
        (28, 20, 0.413556, 482, 108.9, 101.6),
        (29, 21, 0.409635, 487, 107.9, 101.6),
    ];

    fn paper_regime(row_log_size: usize) -> QaRegime {
        QaRegime {
            row_log_size,
            index: 2,
            field_bits: 127,
            distance_failure_bits: 103,
            query_security_bits: 103,
        }
    }

    /// BaseFold query error `(3/4)^249` — the paper's rate-1/2
    /// unique-decoding choice.
    fn paper_basefold_error() -> ErrorBits {
        ErrorBits::from_log2(-249.0 * log2(0.75))
    }

    #[test]
    fn distance_and_queries_match_paper_table_4() {
        for &(_k, n2, delta, t, _, _) in &TABLE_4 {
            let regime = paper_regime(n2);
            let d = regime.distance().expect("Mersenne127 admits a bound");
            assert!(
                libm::fabs(d - delta) < 5e-7,
                "delta mismatch at n2={n2}: got {d}, paper {delta}"
            );
            assert_eq!(regime.num_queries(), Some(t), "query count at n2={n2}");
        }
    }

    #[test]
    fn union_bound_matches_paper_table_4() {
        for &(_k, n2, _, t, lambda_fld, lambda_total) in &TABLE_4 {
            let regime = paper_regime(n2);

            // Each principal term is targeted at ~103 bits.
            assert!(libm::fabs(regime.distance_error().unwrap().bits() - 103.0) < 0.05);
            assert!(libm::fabs(regime.sampling_error(t).unwrap().bits() - 103.1) < 0.15);

            let fld = regime.field_error(t).unwrap().bits();
            assert!(
                libm::fabs(fld - lambda_fld) < 0.05,
                "lambda_fld at n2={n2}: got {fld}, paper {lambda_fld}"
            );

            let total = regime
                .total_error(paper_basefold_error(), ErrorBits::from_log2(128.0))
                .unwrap()
                .bits();
            assert!(
                libm::fabs(total - lambda_total) < 0.05,
                "lambda_total at n2={n2}: got {total}, paper {lambda_total}"
            );
            // The paper's claim: every row clears 100 bits.
            assert!(total >= 100.0);
        }
    }

    /// The auxiliary term is far from binding at these parameters — worth
    /// pinning so a future change that makes it binding is visible.
    #[test]
    fn auxiliary_term_is_slack() {
        let regime = paper_regime(12);
        assert!(regime.auxiliary_error(445).unwrap().bits() > 200.0);
    }

    /// 31-bit fields admit no useful QA distance at any practical length:
    /// `(1 + log₂N/log₂p)/c` already exceeds `g_p(δ) ≤ 1`.
    #[test]
    fn small_fields_have_no_admissible_distance() {
        for row_log_size in [10, 16, 20] {
            for index in [2, 4] {
                let regime = QaRegime {
                    row_log_size,
                    index,
                    field_bits: 31,
                    distance_failure_bits: 100,
                    query_security_bits: 100,
                };
                assert_eq!(regime.distance(), None);
                assert_eq!(regime.num_queries(), None);
                assert_eq!(regime.distance_error(), None);
            }
        }
    }

    /// ed25519 (255-bit) is the operating point for this fork: a better
    /// distance than the paper's Mersenne127, hence fewer opened columns.
    #[test]
    fn ed25519_beats_mersenne127() {
        // Same targets as the paper's Table 4, so only the field differs.
        let make = |field_bits| QaRegime {
            row_log_size: 16,
            index: 2,
            field_bits,
            distance_failure_bits: 103,
            query_security_bits: 103,
        };
        let m127 = make(127);
        let ed = make(255);

        let (d127, ded) = (m127.distance().unwrap(), ed.distance().unwrap());
        assert!(ded > d127, "ed25519 distance {ded} should exceed {d127}");
        assert!(libm::fabs(ded - 0.464720) < 5e-7);
        assert_eq!(ed.num_queries(), Some(425));
        assert!(ed.num_queries().unwrap() < m127.num_queries().unwrap());
    }

    /// Rate 1/4 buys a much larger distance, and so many fewer queries.
    #[test]
    fn higher_index_reduces_queries() {
        let regime = QaRegime {
            row_log_size: 16,
            index: 4,
            field_bits: 255,
            distance_failure_bits: 103,
            query_security_bits: 103,
        };
        assert!(regime.distance().unwrap() > 0.7);
        // vs. 425 at index 2 — a 2× longer codeword for 1.7× fewer columns.
        assert_eq!(regime.num_queries(), Some(256));
    }

    /// `δ` shrinks monotonically as the code gets longer at fixed field
    /// size, so the query count grows monotonically.
    #[test]
    fn queries_grow_with_length() {
        let mut previous = 0;
        for n2 in 12..=21 {
            let t = paper_regime(n2).num_queries().unwrap();
            assert!(t >= previous, "query count should not decrease at n2={n2}");
            previous = t;
        }
    }
}
