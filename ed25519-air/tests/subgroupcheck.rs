//! Structural tests for `EcSubgroupCheckAir` over `p3-bn254`'s field.
//!
//! `p3-bn254` is used here purely as a convenient, already-implemented native field to
//! validate the *shape* of the constraints (helper-column wiring, degree, sub-air composition,
//! etc). The `d` curve parameter below is an arbitrary field element, not the real ed25519
//! constant — that only becomes meaningful once ed25519's actual base field lands in this repo.
//! Likewise, over this placeholder field/curve there's no real "prime-order subgroup" to speak
//! of; these tests only check that the AIR correctly enforces "the row's `(q_x, q_y)` is `8*P`
//! for the row's witness `P`," which is the constraint that becomes a subgroup-membership
//! certificate once the real curve is wired in (see `docs/plans/SubgroupCheck.md`).
//!
//! The Fiat-Shamir transcript below is a minimal insecure test-only challenger (not a real
//! Poseidon2 sponge) built directly from `Bn254` field arithmetic, mirroring `tests/pointcheck.rs`.

use std::marker::PhantomData;

use p3_bn254::Bn254;
use p3_challenger::{CanObserve, CanSample, CanSampleBits, FieldChallenger};
use p3_commit::testing::TrivialPcs;
use p3_dft::Radix2DitParallel;
use num_bigint::BigUint;
use p3_ed25519_air::add::add_points;
use p3_ed25519_air::subgroupcheck::{
    EcSubgroupCheckAir, generate_trace_rows, generate_trace_rows_for_targets, scalar_mul,
};
use p3_field::{Field, PrimeCharacteristicRing, PrimeField};
use p3_uni_stark::{StarkConfig, prove, verify};

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

type Val = Bn254;
type Challenge = Bn254;
type Dft = Radix2DitParallel<Val>;
type Pcs = TrivialPcs<Val, Dft>;
type Config = StarkConfig<Pcs, Challenge, TestChallenger>;

fn build_config(log_n: usize) -> Config {
    let pcs = Pcs {
        dft: Dft::default(),
        log_n,
        _phantom: PhantomData,
    };
    Config::new(pcs, TestChallenger::new())
}

/// Placeholder curve constant; `a = -1` is fixed by `EcSubgroupCheckAir` (matching `EcAddAir`).
/// With `d = 0` the curve is `y^2 - x^2 = 1`.
const D: u64 = 0;

/// On-curve witness points for `a = -1`, `d = 0`, i.e. `y^2 - x^2 = 1`.
///
/// Uses the standard rational parametrization `x = 2t/(1-t^2)`, `y = (1+t^2)/(1-t^2)`, which
/// satisfies `y^2 - x^2 = 1` as a formal identity over any field, as long as `1-t^2` is
/// invertible. `t` starts at 2 (rather than 0) purely to avoid the `t = 1` singularity.
fn sample_points(rows: usize) -> Vec<(Bn254, Bn254)> {
    (0..rows)
        .map(|i| {
            let t = Bn254::from_u32(i as u32 + 2);
            let denom_inv = (Bn254::ONE - t * t).inverse();
            let x = (Bn254::from_u32(2) * t) * denom_inv;
            let y = (Bn254::ONE + t * t) * denom_inv;
            (x, y)
        })
        .collect()
}

#[test]
fn subgroupcheck_prove_and_verify() {
    let d = Bn254::from_u32(D as u32);
    let air = EcSubgroupCheckAir { d };
    let trace = generate_trace_rows(d, &sample_points(8));

    let config = build_config(3);
    let proof = prove(&config, &air, trace, &[]);
    verify(&config, &air, &proof, &[]).expect("verification failed");
}

// In debug builds, `prove` itself checks constraints and panics on a tampered trace,
// before a bad proof could even reach `verify`. See `p3_uni_stark::prove`'s
// `DebugConstraintBuilder` bound and e.g. `uni-stark/tests/fib_air.rs`'s
// `test_incorrect_public_value` for the same pattern.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "constraints not satisfied on row")]
fn subgroupcheck_rejects_mismatched_target() {
    let d = Bn254::from_u32(D as u32);
    let air = EcSubgroupCheckAir { d };
    let mut trace = generate_trace_rows(d, &sample_points(8));

    // Corrupt the first row's `q_x` so it no longer equals `8*P`, leaving the witness `P` and
    // all the doubling arithmetic untouched.
    let width = trace.width;
    let q_x_col = width - 2;
    trace.values[q_x_col] += Bn254::ONE;

    let config = build_config(3);
    let _ = prove(&config, &air, trace, &[]);
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "constraints not satisfied on row")]
fn subgroupcheck_rejects_off_curve_witness() {
    let d = Bn254::from_u32(D as u32);
    let air = EcSubgroupCheckAir { d };
    let mut trace = generate_trace_rows(d, &sample_points(8));

    // Corrupt the first row's witness `p_x` (col 0) and its point-check-block copy `cc_x`
    // (col 2) by the same amount, keeping the `cc_x = p_x` copy constraint satisfied while
    // pushing the witness off-curve, so the embedded `EcPointCheckAir` sub-air's curve
    // equation constraint fails specifically.
    trace.values[0] += Bn254::ONE;
    trace.values[2] += Bn254::ONE;

    let config = build_config(3);
    let _ = prove(&config, &air, trace, &[]);
}

/// `scalar_mul(d, 8, P)` must agree with three explicit doublings via `add_points` — this
/// validates the double-and-add implementation itself, independent of any group-order
/// reasoning (unlike `witness_for_target`, multiplying by the literal integer `8` doesn't
/// depend on knowing the ambient group's order).
#[test]
fn scalar_mul_by_8_matches_tripled_doubling() {
    let d = Bn254::from_u32(D as u32);
    let p = sample_points(1)[0];

    let expected = {
        let p2 = add_points(d, p, p);
        let p4 = add_points(d, p2, p2);
        add_points(d, p4, p4)
    };

    assert_eq!(scalar_mul(d, &BigUint::from(8u32), p), expected);
}

/// `generate_trace_rows_for_targets` derives the witness `P` from a target `Q` internally
/// (rather than requiring the caller to solve for it), given the subgroup's known order `l`.
///
/// This placeholder field/curve has no real prime-order subgroup to pick a genuine `l` from
/// (see the module-level docs), so this only exercises the degenerate but well-defined case
/// `Q = identity`: any `l` divides the identity's order (1), so `8 * witness_for_target(l, O)`
/// always equals `O` again, letting the full derive-witness -> prove -> verify path run
/// end-to-end. A round trip against a non-trivial `Q` needs the real curve's subgroup order,
/// same as the rest of this module's real-vs-placeholder caveat.
#[test]
fn subgroupcheck_for_targets_prove_and_verify_identity() {
    let d = Bn254::from_u32(D as u32);
    let air = EcSubgroupCheckAir { d };
    let identity = (Bn254::ZERO, Bn254::ONE);
    let l = BigUint::from(11u32); // odd prime, coprime to 8; the specific value is irrelevant here
    let trace = generate_trace_rows_for_targets(d, &l, &[identity; 8]);

    let config = build_config(3);
    let proof = prove(&config, &air, trace, &[]);
    verify(&config, &air, &proof, &[]).expect("verification failed");
}
