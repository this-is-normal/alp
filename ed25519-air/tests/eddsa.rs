//! Structural tests for the monolithic `EdDsaAir` over `p3-bn254`'s field.
//!
//! See `tests/subgroupcheck.rs` for why `p3-bn254` and the `d` constant here are placeholders
//! rather than the real ed25519 base field / curve parameter, and why the Fiat-Shamir transcript
//! below is a minimal insecure test-only challenger.
//!
//! **The positive case is deliberately degenerate.** Over this placeholder curve
//! `EcSubgroupCheckAir`'s witness derivation only closes for `target = identity` (see
//! `subgroupcheck_for_targets_prove_and_verify_identity`), so both subgroup-checked points are
//! forced to `A = R = O`, which in turn forces `s = 0` — otherwise `[s]B != R + [k]A` and the
//! verification equation genuinely fails. `k` and the generator `B` stay non-trivial, so every
//! constraint path is still *exercised*; this is a structural test, not a real signature. A
//! non-degenerate positive test needs the real curve (see `docs/EdDSA.md`).
//!
//! The negative tests are the substance here. The circuit's security rests entirely on the
//! shared columns tying the four row types together (`EdDsaSharedCols`), and each tie is a
//! constraint that could be dropped without any *positive* test noticing. So rather than
//! `should_panic`, they use `p3_air::check_all_constraints` and assert **which rows** fail:
//! a test that only checked "something was rejected" would still pass if the tamper tripped some
//! unrelated constraint instead of the tie it is aiming at.

use std::marker::PhantomData;

use num_bigint::BigUint;
use p3_air::{
    AirLayout, check_all_constraints, get_max_constraint_degree, get_symbolic_constraints,
};
use p3_bn254::Bn254;
use p3_challenger::{CanObserve, CanSample, CanSampleBits, FieldChallenger};
use p3_commit::testing::TrivialPcs;
use p3_dft::Radix2DitParallel;
use p3_ed25519_air::add::{Point, add_points};
use p3_ed25519_air::eddsa::{
    EdDsaAir, FIXED_LADDER_START, NUM_EDDSA_COLS, VAR_LADDER_START, add_cols_mut, add_row,
    generate_trace_rows, public_values, sgc_a_row, sgc_cols, sgc_cols_mut, sgc_r_row,
    shared_cols_mut, trace_height,
};
use p3_ed25519_air::subgroupcheck;
use p3_field::{Field, PrimeCharacteristicRing, PrimeField};
use p3_matrix::dense::RowMajorMatrix;
use p3_uni_stark::{
    StarkConfig, prove_with_preprocessed, setup_preprocessed, verify_with_preprocessed,
};

#[derive(Clone, Debug)]
struct TestChallenger {
    state: Bn254,
}

impl TestChallenger {
    const fn new() -> Self {
        Self { state: Bn254::ZERO }
    }
}

impl CanObserve<Bn254> for TestChallenger {
    fn observe(&mut self, value: Bn254) {
        self.state = self.state * Bn254::from_u32(1_000_003) + value + Bn254::ONE;
    }
}

impl CanObserve<Vec<Vec<Bn254>>> for TestChallenger {
    fn observe(&mut self, value: Vec<Vec<Bn254>>) {
        for row in value {
            for v in row {
                CanObserve::<Bn254>::observe(self, v);
            }
        }
    }
}

impl CanSample<Bn254> for TestChallenger {
    fn sample(&mut self) -> Bn254 {
        self.state = self.state * Bn254::from_u32(2_000_029) + Bn254::ONE;
        self.state
    }
}

impl CanSampleBits<usize> for TestChallenger {
    fn sample_bits(&mut self, bits: usize) -> usize {
        assert!(bits < usize::BITS as usize);
        let sampled = CanSample::<Bn254>::sample(self);
        let bytes = sampled.as_canonical_biguint().to_bytes_le();
        let mut buf = [0u8; 8];
        let n = bytes.len().min(8);
        buf[..n].copy_from_slice(&bytes[..n]);
        let val = u64::from_le_bytes(buf) as usize;
        let mask = (1usize << bits) - 1;
        val & mask
    }
}

impl FieldChallenger<Bn254> for TestChallenger {}

type F = Bn254;
type Challenge = Bn254;
type Dft = Radix2DitParallel<F>;
type Pcs = TrivialPcs<F, Dft>;
type Config = StarkConfig<Pcs, Challenge, TestChallenger>;

fn build_config(log_n: usize) -> Config {
    let pcs = Pcs {
        dft: Dft::default(),
        log_n,
        _phantom: PhantomData,
    };
    Config::new(pcs, TestChallenger::new())
}

/// Placeholder curve constant; `a = -1` is fixed by the AIR. With `d = 0` the curve is
/// `y^2 - x^2 = 1`, whose rational parametrization gives cheap on-curve points ([`curve_point`]).
const D: u64 = 0;

/// Scalar bit length. Real ed25519 uses 253; five keeps the trace at eight rows.
const N: usize = 5;

/// `log2(trace_height(N))` = `log2(8)`.
const LOG_HEIGHT: usize = 3;

fn d() -> F {
    F::from_u64(D)
}

/// The neutral element, on-curve for every `d`.
const fn identity() -> Point<F> {
    (F::ZERO, F::ONE)
}

/// On-curve points for `a = -1`, `d = 0`, i.e. `y^2 - x^2 = 1`, via the standard rational
/// parametrization `x = 2t/(1-t^2)`, `y = (1+t^2)/(1-t^2)`. `t` starts at 2 to avoid the `t = 1`
/// singularity. Same construction as `tests/subgroupcheck.rs`.
fn curve_point(i: u32) -> Point<F> {
    let t = F::from_u32(i + 2);
    let denom_inv = (F::ONE - t * t).inverse();
    ((F::TWO * t) * denom_inv, (F::ONE + t * t) * denom_inv)
}

/// The generator `B`. Unlike `A` and `R` it is never subgroup-checked — it is baked into the
/// preprocessed trace — so it can be a non-trivial point even in the degenerate setting.
fn base_b() -> Point<F> {
    curve_point(0)
}

/// Subgroup order, used only to derive the subgroup-check witnesses. Any odd value coprime to 8
/// works for `target = identity`; see `tests/subgroupcheck.rs`.
fn l() -> BigUint {
    BigUint::from(11u32)
}

/// LSB-first bit decomposition of `scalar` into exactly [`N`] bits.
fn bits(scalar: u32) -> Vec<bool> {
    assert!(scalar < 1 << N);
    (0..N).map(|i| (scalar >> i) & 1 == 1).collect()
}

/// One (public key, signature, challenge) instance, and the trace/AIR/public values it induces.
#[derive(Clone)]
struct Scenario {
    a: Point<F>,
    r: Point<F>,
    s_bits: Vec<bool>,
    k_bits: Vec<bool>,
}

impl Scenario {
    /// The satisfying instance: `A = R = O` forces `s = 0`, leaving `k` free. See the module docs
    /// for why it cannot be less degenerate over the placeholder curve.
    fn valid() -> Self {
        Self {
            a: identity(),
            r: identity(),
            s_bits: bits(0),
            k_bits: bits(0b0_1101),
        }
    }

    fn air(&self) -> EdDsaAir<F> {
        EdDsaAir {
            d: d(),
            base_b: base_b(),
            num_bits: N,
        }
    }

    fn trace(&self) -> RowMajorMatrix<F> {
        generate_trace_rows(
            d(),
            base_b(),
            &l(),
            self.a,
            self.r,
            &self.s_bits,
            &self.k_bits,
        )
    }

    fn public_values(&self) -> Vec<F> {
        public_values(self.a, self.r, &self.s_bits, &self.k_bits)
    }
}

/// The distinct trace rows on which at least one constraint fails, in ascending order.
///
/// Asserting on rows (rather than just "it was rejected") is what makes the negative tests below
/// specific: each tamper names the row whose tie it targets.
fn failing_rows(air: &EdDsaAir<F>, trace: &RowMajorMatrix<F>, public_values: &[F]) -> Vec<usize> {
    let report = check_all_constraints(air, trace, public_values, None);
    let mut rows: Vec<usize> = report.failures.iter().map(|f| f.row).collect();
    // `check_all_constraints` walks rows in order, so equal rows are already adjacent.
    rows.dedup();
    rows
}

/// Mutable view of one row of a main trace.
fn row_mut(trace: &mut RowMajorMatrix<F>, row: usize) -> &mut [F] {
    &mut trace.values[row * NUM_EDDSA_COLS..(row + 1) * NUM_EDDSA_COLS]
}

// -- positive ------------------------------------------------------------------------------------

#[test]
fn eddsa_prove_and_verify() {
    let scenario = Scenario::valid();
    let air = scenario.air();
    let trace = scenario.trace();
    let public_values = scenario.public_values();

    let config = build_config(LOG_HEIGHT);
    let (preprocessed_prover_data, preprocessed_vk) =
        setup_preprocessed::<Config, _>(&config, &air, LOG_HEIGHT).unwrap();

    let proof = prove_with_preprocessed(
        &config,
        &air,
        trace,
        &public_values,
        Some(&preprocessed_prover_data),
    );
    verify_with_preprocessed(
        &config,
        &air,
        &proof,
        &public_values,
        Some(&preprocessed_vk),
    )
    .expect("verification failed");
}

/// The same instance, checked constraint-by-constraint rather than through the proof system: no
/// row may fail, including the aux rows and any padding.
#[test]
fn eddsa_valid_trace_satisfies_every_constraint() {
    let scenario = Scenario::valid();
    let trace = scenario.trace();
    assert_eq!(trace.values.len() / NUM_EDDSA_COLS, trace_height(N));
    assert!(failing_rows(&scenario.air(), &trace, &scenario.public_values()).is_empty());
}

/// Guards the cost claim documented in `docs/EdDSA.md`. Degree 3 is the load-bearing half: the
/// building blocks are degree 2 and one selector factor makes them 3, so a regression here means
/// a constraint picked up an extra witness factor and the proof got more expensive.
#[test]
fn eddsa_constraint_count_and_degree() {
    let air = Scenario::valid().air();
    let layout = AirLayout::from_air(&air);
    assert_eq!(get_symbolic_constraints::<F, _>(&air, layout).len(), 116);
    assert_eq!(
        get_max_constraint_degree::<F, _>(&air, layout, trace_height(N)),
        3
    );
}

// -- negative: one per tie between row types -----------------------------------------------------
//
// Every test below starts from the satisfying instance and breaks exactly one link, then asserts
// the failure lands on the row that owns that link.

/// The subgroup-check row for `A` proves `8P = Q` for its own witness `P` — but nothing in that
/// block says `Q` is the `A` the ladder uses. Here the row is replaced with an *internally valid*
/// subgroup check on a different target, so the only thing left to catch it is
/// `is_sgc_a => (q_x, q_y) == (a_x, a_y)`.
#[test]
fn eddsa_rejects_subgroup_check_of_wrong_target() {
    let scenario = Scenario::valid();
    let mut trace = scenario.trace();

    // A self-consistent subgroup-check row for the witness `curve_point(3)`, whose target `8P`
    // is not the identity and therefore not `A`.
    let other = subgroupcheck::generate_trace_rows(d(), &[curve_point(3)]);
    *sgc_cols_mut(row_mut(&mut trace, sgc_a_row(N))) = *sgc_cols(&other.values);

    assert_eq!(
        failing_rows(&scenario.air(), &trace, &scenario.public_values()),
        vec![sgc_a_row(N)]
    );
}

/// `EcScalarMulAir` treats its base as a pure witness, so on its own the variable ladder would
/// happily compute `[k]P` for any `P` of the prover's choosing — including one outside the
/// subgroup. Row 0's `(x_dbl, y_dbl) == (a_x, a_y)` is what closes that gap.
///
/// Tampering a single column would not test that: the row-0 add law derives `sum` from
/// `(x_dbl, y_dbl)`, so a lone edit trips the law and the test would pass even with the tie
/// deleted. Instead the *entire* variable ladder is replaced with a self-consistent one over a
/// different base, chosen so that even the downstream `[k]A` readout still matches — leaving the
/// row-0 tie as the only constraint standing.
#[test]
fn eddsa_rejects_variable_ladder_base_not_a() {
    // `(0, -1)` has order 2 (see the assertion below), so for an *even* `k` the substituted
    // ladder still ends at `[k]P = O = [k]A` and the `is_ladder_last` readout tie is satisfied.
    let p = (F::ZERO, F::NEG_ONE);
    assert_eq!(add_points(d(), p, p), identity(), "P must have order 2");

    let scenario = Scenario {
        k_bits: bits(0b0_1100), // 12, even
        ..Scenario::valid()
    };
    let mut trace = scenario.trace();

    // Identical instance except for the ladder's base, so the substituted rows agree with the
    // host trace on the scalar bits, the weights and the accumulator's final value.
    let substitute = Scenario {
        a: p,
        ..scenario.clone()
    }
    .trace();
    for row in 0..N {
        let lo = row * NUM_EDDSA_COLS + VAR_LADDER_START;
        let hi = row * NUM_EDDSA_COLS + FIXED_LADDER_START;
        trace.values[lo..hi].copy_from_slice(&substitute.values[lo..hi]);
    }

    assert_eq!(
        failing_rows(&scenario.air(), &trace, &scenario.public_values()),
        vec![0]
    );
}

/// `[k]A` is a witness, not a public value, so a prover is free to claim a different one — as long
/// as the claim is consistent everywhere it is read. Rewriting the shared `ka_x` on *every* row
/// keeps the propagation constraint happy, and leaves exactly the two rows that consume it: the
/// ladder's last row, which must produce it, and the add row, which must add it to `R`.
#[test]
fn eddsa_rejects_claimed_ka_not_matching_ladder_result() {
    let scenario = Scenario::valid();
    let mut trace = scenario.trace();

    for row in 0..trace_height(N) {
        shared_cols_mut(row_mut(&mut trace, row)).ka_x += F::ONE;
    }

    assert_eq!(
        failing_rows(&scenario.air(), &trace, &scenario.public_values()),
        vec![N - 1, add_row(N)]
    );
}

/// The add row's first operand must be the subgroup-checked `R`, not an arbitrary point.
///
/// Unlike the tests around it, this one pins the *row* rather than the individual tie, and that
/// is not a gap in the test — the tie cannot be violated alone. Edwards addition is a bijection
/// in each argument, so once `(x2, y2) == ka` and `(x3, y3) == sb` hold, the add law determines
/// `(x1, y1)` outright; any trace with a different first operand also breaks the law or the final
/// equality. (Confirmed by mutation: deleting `x1 == r_x` alone still leaves this rejected.) What
/// the tie buys is not extra rejection on this row but *direction* — it forces the determined
/// operand to be the `R` that row `sgc_r` checked, rather than leaving `R` unconstrained.
#[test]
fn eddsa_rejects_add_row_first_operand_not_r() {
    let scenario = Scenario::valid();
    let mut trace = scenario.trace();

    add_cols_mut(row_mut(&mut trace, add_row(N))).x1 += F::ONE;

    assert_eq!(
        failing_rows(&scenario.air(), &trace, &scenario.public_values()),
        vec![add_row(N)]
    );
}

/// The forgery the whole circuit exists to reject: a *well-formed* trace for an instance where
/// the verification equation simply does not hold. `s = 1` against `A = R = O` gives
/// `[s]B = B != O = R + [k]A`. Nothing is tampered — both ladders, both subgroup checks and the
/// addition are all internally correct — so only the add row's final `(x3, y3) == (sb_x, sb_y)`
/// stands between this and a valid proof.
#[test]
fn eddsa_rejects_signature_failing_the_verification_equation() {
    let forged = Scenario {
        s_bits: bits(1),
        ..Scenario::valid()
    };
    let trace = forged.trace();

    assert_eq!(
        failing_rows(&forged.air(), &trace, &forged.public_values()),
        vec![add_row(N)]
    );
}

/// The shared columns are only a valid wire if they cannot change mid-table: otherwise a prover
/// could satisfy the ladder's readout with one `s` and the first-row public binding with another.
/// `s` is read on row 0 and row `N - 1` only, so breaking propagation at the `sgc_a -> sgc_r`
/// transition — where no other constraint reads it — isolates the propagation constraint itself.
#[test]
fn eddsa_rejects_shared_column_changing_mid_trace() {
    let scenario = Scenario::valid();
    let mut trace = scenario.trace();

    for row in sgc_r_row(N)..trace_height(N) {
        shared_cols_mut(row_mut(&mut trace, row)).s += F::ONE;
    }

    assert_eq!(
        failing_rows(&scenario.air(), &trace, &scenario.public_values()),
        vec![sgc_a_row(N)]
    );
}

/// The public values are bound to the shared columns on the first row and nowhere else; that
/// single binding, plus propagation, is what makes the whole trace answer to the *claimed*
/// `(A, R, s, k)` rather than to whatever the prover found convenient.
#[test]
fn eddsa_rejects_public_value_not_matching_row_zero() {
    let scenario = Scenario::valid();
    let trace = scenario.trace();

    let mut public_values = scenario.public_values();
    public_values[4] += F::ONE; // `s`

    assert_eq!(
        failing_rows(&scenario.air(), &trace, &public_values),
        vec![0]
    );
}

/// The proof system itself must reject a tampered trace, not just the standalone constraint
/// checker. In debug builds `prove_with_preprocessed` checks constraints internally and panics
/// before a bad proof could reach the verifier; see `tests/subgroupcheck.rs` for the same pattern.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "constraints not satisfied on row")]
fn eddsa_prover_rejects_tampered_trace() {
    let scenario = Scenario::valid();
    let air = scenario.air();
    let mut trace = scenario.trace();
    add_cols_mut(row_mut(&mut trace, add_row(N))).x3 += F::ONE;

    let config = build_config(LOG_HEIGHT);
    let (preprocessed_prover_data, _) =
        setup_preprocessed::<Config, _>(&config, &air, LOG_HEIGHT).unwrap();
    let _ = prove_with_preprocessed(
        &config,
        &air,
        trace,
        &scenario.public_values(),
        Some(&preprocessed_prover_data),
    );
}
