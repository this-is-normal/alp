//! Structural tests for `EcPointCheckAir` over `p3-bn254`'s field.
//!
//! `p3-bn254` is used here purely as a convenient, already-implemented native field to
//! validate the *shape* of the curve-membership constraints (helper-column wiring, degree,
//! etc). `a = -1` is fixed (matching ed25519); the `d` curve parameter below is an arbitrary
//! field element, not the real ed25519 constant — that only becomes meaningful once
//! ed25519's actual base field lands in this repo.
//!
//! The Fiat-Shamir transcript below is a minimal insecure test-only challenger (not a real
//! Poseidon2 sponge) built directly from `Bn254` field arithmetic, since no ready-made
//! challenger exists for this field in the repo yet.

use std::marker::PhantomData;

use p3_bn254::Bn254;
use p3_challenger::{CanObserve, CanSample, CanSampleBits, FieldChallenger};
use p3_commit::testing::TrivialPcs;
use p3_dft::Radix2DitParallel;
use p3_ed25519_air::pointcheck::{EcPointCheckAir, generate_trace_rows};
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

/// Placeholder curve constant; see module docs. `a = -1` is fixed (matching ed25519 and
/// `EcAddAir`). Chosen so `(0, 1)` (the padding identity) and the sample points below all
/// satisfy `-x^2 + y^2 = 1 + d*x^2*y^2`.
const D: u64 = 0;

/// On-curve points for `a = -1`, `d = 0`, i.e. `y^2 - x^2 = 1`.
///
/// Uses the standard rational parametrization of the unit hyperbola,
/// `x = 2t/(1-t^2)`, `y = (1+t^2)/(1-t^2)`, which satisfies `y^2 - x^2 = 1` as a formal
/// identity over any field (no square roots needed), as long as `1-t^2` is invertible.
/// `t` starts at 2 (rather than 0) so no sample hits `t^2 = 1`, where the denominator
/// vanishes.
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
fn pointcheck_prove_and_verify() {
    let d = Bn254::from_u32(D as u32);
    let air = EcPointCheckAir { d };
    let trace = generate_trace_rows(&sample_points(8));

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
fn pointcheck_rejects_off_curve_point() {
    let d = Bn254::from_u32(D as u32);
    let air = EcPointCheckAir { d };
    let mut trace = generate_trace_rows(&sample_points(8));

    // Corrupt the `x` coordinate of the first row so it no longer matches its `x2` helper
    // column (and is off-curve).
    trace.values[0] += Bn254::ONE;

    let config = build_config(3);
    let _ = prove(&config, &air, trace, &[]);
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "constraints not satisfied on row")]
fn pointcheck_rejects_tampered_helper_column() {
    let d = Bn254::from_u32(D as u32);
    let air = EcPointCheckAir { d };
    let mut trace = generate_trace_rows(&sample_points(8));

    // Corrupt the `x2` helper column of the first row, leaving `x`, `y` valid: the
    // `x2 == x*x` consistency constraint should fail even though `x`, `y` are still on-curve.
    trace.values[2] += Bn254::ONE;

    let config = build_config(3);
    let _ = prove(&config, &air, trace, &[]);
}
