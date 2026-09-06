//! Distance-derived query counts and vetted parameter presets for the Brakedown
//! (recursive Spielman expander) code — the real, security-derived replacement for
//! the old `DEMO_QUERIES = 32` constant.
//!
//! # Query count
//! For the tensor/Ligero proximity test at `security_bits`-bit soundness the number
//! of opened columns is
//! ```text
//!   t = ceil( security_bits / -log2(1 - e) )
//! ```
//! where `e` is the tested proximity radius, expressed as a fraction of the code's
//! relative minimum distance `delta`:
//! * **provable** (Ligero / Roth–Zémor, general linear codes): `e = delta/3`
//! * **conjectured** (Diamond–Posen *"Conjecture 1"*, the Reed–Solomon-proven
//!   analogue up to the unique-decoding radius): `e = delta/2`
//!
//! `delta = beta / r` is the code's (assumed) relative minimum distance. References:
//! Haböck 2023 *"Brakedown's expander code"* (the `(alpha, rho, delta)` sets and the
//! degree derivation), and `conroi/lcpc` (whose reference impl uses `e = delta/3`).
//!
//! # Presets
//! The three presets are the standard operating points from Haböck's Table 1
//! (= lcpc `SdigCode1/4/6`), row weights `c_a = c_n`, `c_b = d_n`. They trade proof
//! size against prover (encode) cost along the code's distance axis:
//! * [`small_proof`]  `delta=0.070` — fewest queries (smallest proof), heaviest encode.
//! * [`balanced`]     `delta=0.050` — the middle.
//! * [`fast_prover`]  `delta=0.020` — cheapest encode, most queries (largest proof).
//!
//! The row weights are Haböck's 128-bit / λ=110 values; over a 255-bit field at
//! λ=100 the required degrees are *lower*, so these are safe over-estimates.

use libm::log2;

use crate::BrakedownParams;

/// A named Brakedown operating point: code parameters plus the (assumed) relative
/// minimum distance `delta = beta/r` used to derive the query count.
#[derive(Clone, Debug)]
pub struct Preset {
    pub name: &'static str,
    /// Assumed relative minimum distance `delta = beta/r` (Haböck Table 1 / lcpc).
    pub delta: f64,
    pub params: BrakedownParams,
}

impl Preset {
    /// Column-query count for `security_bits`-bit soundness in the chosen regime
    /// (`conjectured = false` → provable `delta/3`; `true` → conjectured `delta/2`).
    #[must_use]
    pub fn num_queries(&self, security_bits: usize, conjectured: bool) -> usize {
        brakedown_num_queries(self.delta, security_bits, conjectured)
    }
}

/// `t = ceil( security_bits / -log2(1 - e) )`, with `e = delta/3` (provable) or
/// `e = delta/2` (conjectured). `no_std` ceil via integer comparison.
#[must_use]
pub fn brakedown_num_queries(delta: f64, security_bits: usize, conjectured: bool) -> usize {
    let e = if conjectured { delta / 2.0 } else { delta / 3.0 };
    let per_query = -log2(1.0 - e); // bits of soundness detected per opened column
    let t = security_bits as f64 / per_query;
    let floor = t as usize;
    if (floor as f64) < t { floor + 1 } else { floor }
}

/// Small-proof point — Haböck `SdigCode6`: `alpha=0.238`, `rho=0.581`, `delta=0.070`,
/// row weights `9/12`. Fewest queries, heaviest linear encode.
#[must_use]
pub fn small_proof() -> Preset {
    Preset {
        name: "small-proof",
        delta: 0.070,
        params: BrakedownParams {
            base_threshold: 16,
            alpha_num: 238,
            alpha_den: 1000,
            rate_num: 172,
            rate_den: 100,
            c_a: 9,
            c_b: 12,
        },
    }
}

/// Balanced point — Haböck `SdigCode4`: `alpha=0.200`, `rho=0.610`, `delta=0.050`,
/// row weights `7/12`.
#[must_use]
pub fn balanced() -> Preset {
    Preset {
        name: "balanced",
        delta: 0.050,
        params: BrakedownParams {
            base_threshold: 16,
            alpha_num: 200,
            alpha_den: 1000,
            rate_num: 164,
            rate_den: 100,
            c_a: 7,
            c_b: 12,
        },
    }
}

/// Fast-prover point — Haböck `SdigCode1`: `alpha=0.120`, `rho=0.704`, `delta=0.020`,
/// row weights `6/16`. Cheapest encode, most queries.
#[must_use]
pub fn fast_prover() -> Preset {
    Preset {
        name: "fast-prover",
        delta: 0.020,
        params: BrakedownParams {
            base_threshold: 16,
            alpha_num: 120,
            alpha_den: 1000,
            rate_num: 142,
            rate_den: 100,
            c_a: 6,
            c_b: 16,
        },
    }
}

/// The three presets in order: small-proof → balanced → fast-prover.
#[must_use]
pub fn all() -> [Preset; 3] {
    [small_proof(), balanced(), fast_prover()]
}

/// Optimal row/column split for the tensor commitment (the paper's dimension choice, vs the naive
/// square `log_rows = num_vars/2`). Chooses `log_rows` to minimise the proof term `t·rows + cols`
/// (one combined row, after the Diamond–Posen test/eval consolidation) over the power-of-2 aspect
/// ratios whose codeword has at least `t` columns — so `t` queries (hence the target security)
/// actually fit. Wider matrices ⇒ shorter opened columns ⇒ far smaller proof and verify.
#[must_use]
pub fn optimal_log_rows(num_vars: usize, t: usize, params: &BrakedownParams) -> usize {
    let mut best_lr = 0usize;
    let mut best_cost = u128::MAX;
    for lr in 0..=num_vars {
        let cols = 1u128 << (num_vars - lr);
        // codeword_len ≈ rate · cols  (rate = rate_num/rate_den, the inverse rate).
        let codeword = cols * params.rate_num as u128 / params.rate_den as u128;
        if codeword < t as u128 {
            continue; // too few columns to run t queries → can't reach target security
        }
        let cost = t as u128 * (1u128 << lr) + cols; // proof ∝ t·rows + cols
        if cost < best_cost {
            best_cost = cost;
            best_lr = lr;
        }
    }
    best_lr
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Query counts match `t = ceil(100 / -log2(1 - e))` for both regimes.
    /// (±1 tolerates float-ceil boundary jitter across platforms.)
    #[test]
    fn query_counts_match_derivation() {
        let approx = |got: usize, want: i64| {
            assert!((got as i64 - want).abs() <= 1, "got {got}, want ~{want}");
        };
        // small-proof  delta=0.070
        approx(brakedown_num_queries(0.070, 100, false), 2936); // provable  delta/3
        approx(brakedown_num_queries(0.070, 100, true), 1946); //  conjectured delta/2
        // balanced     delta=0.050
        approx(brakedown_num_queries(0.050, 100, false), 4125);
        approx(brakedown_num_queries(0.050, 100, true), 2738);
        // fast-prover  delta=0.020
        approx(brakedown_num_queries(0.020, 100, false), 10363);
        approx(brakedown_num_queries(0.020, 100, true), 6897);
    }

    /// The distance axis is a genuine tradeoff, and conjectured is always cheaper.
    #[test]
    fn monotonic_tradeoff() {
        let [s, b, f] = all();
        // higher delta ⇒ fewer queries: small-proof < balanced < fast-prover
        assert!(s.num_queries(100, false) < b.num_queries(100, false));
        assert!(b.num_queries(100, false) < f.num_queries(100, false));
        // conjectured (delta/2) always fewer than provable (delta/3)
        for p in [&s, &b, &f] {
            assert!(p.num_queries(100, true) < p.num_queries(100, false));
        }
        // distances and degrees are as specified
        assert!(s.delta > b.delta && b.delta > f.delta);
        assert_eq!((s.params.c_a, s.params.c_b), (9, 12));
        assert_eq!((f.params.c_a, f.params.c_b), (6, 16));
        // sanity: real security is far above the old demo constant of 32
        assert!(s.num_queries(100, true) > 1000);
    }
}
