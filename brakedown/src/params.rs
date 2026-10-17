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

use libm::{ceil, log2, round};

use crate::BrakedownParams;

/// Serde framing beyond the raw element count; calibrated to the consolidated+optimal impl
/// (mirrors `PROOF_OVERHEAD` in `brakedown_cost_model.py`). Only affects `MinSize` scoring.
const PROOF_OVERHEAD: f64 = 1.14;

/// A named Brakedown operating point: code parameters plus the (assumed) relative
/// minimum distance `delta = beta/r` used to derive the query count.
#[derive(Clone, Debug)]
pub struct Preset {
    pub name: &'static str,
    /// Assumed relative minimum distance `delta = beta/r` (Haböck Table 1 / lcpc).
    pub delta: f64,
    /// Linear-encode field-ops per element `Tₙ/n` (verifier re-encode cost; Haböck Table 1).
    pub tn_over_n: f64,
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
        tn_over_n: 18.3,
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
        tn_over_n: 13.7,
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
        tn_over_n: 9.9,
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

/// Row/column split strategy for the tensor commitment. Mirror of the `split` options in
/// `brakedown_cost_model.py`. All strategies keep `codeword_len ≥ t` (enough columns to run `t`
/// queries → the target security is reachable), replacing the naive square split.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SplitStrategy {
    /// Square-ish `log_rows = num_vars/2` (the old default; insecure at small N).
    SquareMatrix,
    /// Closed-form balanced-proof rule `rows ≈ sqrt(N/t)` (Ligero/Brakedown).
    Paper,
    /// Minimise the actual proof bytes.
    MinSize,
    /// Minimise verifier field ops (MUL + ADD).
    MinArith,
    /// Minimise verifier hash (Poseidon2 perm) count.
    MinHashes,
    /// Minimise weighted guest verify cost `hash_cost·hashes + op_cost·arith`.
    MinVerify,
}

/// Hash/serialization cost model used to score a split (mirror of the Python `HASHES` + weights).
#[derive(Clone, Copy, Debug)]
pub struct SplitModel {
    /// Poseidon2 perms to absorb one committed field element (the hash rate).
    pub perms_per_elem: f64,
    /// Merkle digest size, in hash-field elements.
    pub digest_elems: usize,
    /// Bytes per hash-field element.
    pub field_bytes: usize,
    /// Bytes per committed field element (32 for ed25519).
    pub data_elem_bytes: usize,
    /// `MinVerify` weights: guest cost per verifier hash and per verifier field op.
    pub hash_cost: f64,
    pub op_cost: f64,
}

impl SplitModel {
    /// ed25519 data + BabyBear-Poseidon2 (`bb-hash`) Merkle, equal `MinVerify` weights.
    #[must_use]
    pub fn bb_hash() -> Self {
        Self { perms_per_elem: 9.0 / 8.0, digest_elems: 8, field_bytes: 4,
               data_elem_bytes: 32, hash_cost: 1.0, op_cost: 1.0 }
    }
    /// ed25519 data + Poseidon2-over-ed25519 (native) Merkle, equal `MinVerify` weights.
    #[must_use]
    pub fn native() -> Self {
        Self { perms_per_elem: 0.5, digest_elems: 1, field_bytes: 32,
               data_elem_bytes: 32, hash_cost: 1.0, op_cost: 1.0 }
    }
}

/// `(proof_bytes, verifier_hashes, verifier_arith)` for split `log_rows = lr` at query count `t`.
/// COMBOS = 1 (Diamond–Posen consolidation). Matches `split_metrics` in the Python model.
fn split_costs(num_vars: usize, lr: usize, t: usize, preset: &Preset, model: &SplitModel) -> (f64, f64, f64) {
    let rows = (1u128 << lr) as f64;
    let cols = (1u128 << (num_vars - lr)) as f64;
    let codeword = cols * preset.params.rate_num as f64 / preset.params.rate_den as f64;
    let mut depth = ceil(log2(codeword));
    if depth < 1.0 {
        depth = 1.0;
    }
    let leaf = ceil(rows * model.perms_per_elem); // per-query leaf hash (a rows-tall column)
    let tf = t as f64;
    let field_elems = tf * rows + cols;
    let hash_codes = tf * depth;
    let proof_bytes = (field_elems * model.data_elem_bytes as f64
        + hash_codes * model.digest_elems as f64 * model.field_bytes as f64)
        * PROOF_OVERHEAD;
    let v_hashes = tf * (leaf + depth) + cols * model.perms_per_elem; // leaves+paths + FS absorb
    let v_arith = 2.0 * (tf * rows + (preset.tn_over_n / 2.0 + 1.0) * cols); // MUL + ADD
    (proof_bytes, v_hashes, v_arith)
}

/// Choose `log_rows` for the tensor split under `strategy` (replaces the old `optimal_log_rows`,
/// which was `MinSize`). See [`SplitStrategy`]. `MinSize`/`MinArith`/`MinHashes`/`MinVerify`
/// brute-force the power-of-2 aspect ratios; `SquareMatrix`/`Paper` are direct formulas.
#[must_use]
pub fn choose_log_rows(
    num_vars: usize,
    t: usize,
    preset: &Preset,
    model: &SplitModel,
    strategy: SplitStrategy,
) -> usize {
    let enough = |lr: usize| -> bool {
        let cols = (1u128 << (num_vars - lr)) as f64;
        cols * preset.params.rate_num as f64 / preset.params.rate_den as f64 >= t as f64
    };

    match strategy {
        SplitStrategy::SquareMatrix => return num_vars / 2,
        SplitStrategy::Paper => {
            // rows ≈ sqrt(N/t): log_rows = round((num_vars - log2 t)/2), widened until enough.
            let raw = round((num_vars as f64 - log2(t.max(1) as f64)) / 2.0);
            let mut lr = (raw.max(0.0) as usize).min(num_vars);
            while lr > 0 && !enough(lr) {
                lr -= 1;
            }
            assert!(
                enough(lr),
                "Brakedown Paper: no split of a 2^{num_vars} poly can run {t} queries at this \
                 rate — increase num_vars or choose a higher-\u{03b4} preset"
            );
            return lr;
        }
        _ => {}
    }

    let mut best_lr: Option<usize> = None;
    let mut best = f64::INFINITY;
    for lr in 0..=num_vars {
        if !enough(lr) {
            continue; // this shape has too few columns for t queries; try a wider one
        }
        let (pb, vh, va) = split_costs(num_vars, lr, t, preset, model);
        let score = match strategy {
            SplitStrategy::MinSize => pb,
            SplitStrategy::MinArith => va,
            SplitStrategy::MinHashes => vh,
            SplitStrategy::MinVerify => model.hash_cost * vh + model.op_cost * va,
            SplitStrategy::SquareMatrix | SplitStrategy::Paper => unreachable!(),
        };
        if best_lr.is_none() || score < best {
            best = score;
            best_lr = Some(lr);
        }
    }
    // No aspect ratio has enough columns to run t queries → the target security is
    // unreachable at this N; fail loudly rather than silently returning an insecure split.
    best_lr.unwrap_or_else(|| {
        panic!(
            "Brakedown {strategy:?}: no split of a 2^{num_vars} poly can run {t} queries \
             (rate {}/{}) — even the widest codeword has < {t} columns; increase num_vars \
             or choose a higher-\u{03b4} preset",
            preset.params.rate_num, preset.params.rate_den
        )
    })
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

    /// Split strategies land on the Python model's `log_rows` (small-proof, t=2936).
    #[test]
    fn split_strategies() {
        let p = small_proof();
        let m = SplitModel::bb_hash();
        let t = p.num_queries(100, false); // 2936
        // nv=16: square → 256×256 (lr=8); min_size → 4×16384 (lr=2); min_arith is taller.
        assert_eq!(choose_log_rows(16, t, &p, &m, SplitStrategy::SquareMatrix), 8);
        assert_eq!(choose_log_rows(16, t, &p, &m, SplitStrategy::MinSize), 2);
        assert_eq!(choose_log_rows(16, t, &p, &m, SplitStrategy::Paper), 2);
        assert!(choose_log_rows(16, t, &p, &m, SplitStrategy::MinArith) > 2);
        // nv=24: min_size → 64×262144 (lr=6).
        assert_eq!(choose_log_rows(24, t, &p, &m, SplitStrategy::MinSize), 6);
        // MinVerify weighted to arith-only reproduces MinArith; hash-only reproduces MinHashes.
        let arith_only = SplitModel { hash_cost: 0.0, op_cost: 1.0, ..m };
        let hash_only = SplitModel { hash_cost: 1.0, op_cost: 0.0, ..m };
        assert_eq!(
            choose_log_rows(24, t, &p, &arith_only, SplitStrategy::MinVerify),
            choose_log_rows(24, t, &p, &m, SplitStrategy::MinArith)
        );
        assert_eq!(
            choose_log_rows(24, t, &p, &hash_only, SplitStrategy::MinVerify),
            choose_log_rows(24, t, &p, &m, SplitStrategy::MinHashes)
        );
    }

    /// A too-small N where no split can run the required queries must fail loudly, not
    /// silently return an insecure (over-clamped) split.
    #[test]
    #[should_panic(expected = "no split")]
    fn split_unreachable_panics() {
        let p = small_proof();
        let m = SplitModel::bb_hash();
        let t = p.num_queries(100, false); // 2936 — needs >~1700 columns
        // 2^8 = 256 message columns max ⇒ widest codeword ~440 < 2936: impossible.
        let _ = choose_log_rows(8, t, &p, &m, SplitStrategy::MinSize);
    }
}
