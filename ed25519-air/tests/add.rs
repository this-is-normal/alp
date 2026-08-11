//! Structural tests for `EcAddAir` over `p3-bn254`'s field.
//!
//! `p3-bn254` is used here purely as a convenient, already-implemented native field to
//! validate the *shape* of the addition constraints (helper-column wiring, degree, etc).
//! The `d` curve parameter below is an arbitrary field element, not the real ed25519 `d` —
//! that only becomes meaningful once ed25519's actual base field lands in this repo. Since
//! the AIR only checks that `(x3, y3)` is consistent with `(x1, y1)`, `(x2, y2)` under the
//! twisted-Edwards addition law for whatever `d` is supplied, this is still a faithful test
//! of the constraint wiring itself.
//!
//! The Fiat-Shamir transcript below is a minimal insecure test-only challenger (not a real
//! Poseidon2 sponge) built directly from `Bn254` field arithmetic, since no ready-made
//! challenger exists for this field in the repo yet.

use std::marker::PhantomData;

use p3_bn254::Bn254;
use p3_challenger::{CanObserve, CanSample, CanSampleBits, FieldChallenger};
use p3_commit::testing::TrivialPcs;
use p3_dft::Radix2DitParallel;
use p3_ed25519_air::add::{EcAddAir, generate_trace_rows};
use p3_field::{PrimeCharacteristicRing, PrimeField};
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

/// Placeholder curve constant; see module docs.
const D: u64 = 7;

fn sample_pairs(rows: usize) -> Vec<((Bn254, Bn254), (Bn254, Bn254))> {
    (0..rows)
        .map(|i| {
            let i = i as u32;
            (
                (Bn254::from_u32(2 * i + 1), Bn254::from_u32(3 * i + 5)),
                (Bn254::from_u32(4 * i + 2), Bn254::from_u32(5 * i + 11)),
            )
        })
        .collect()
}

#[test]
fn add_prove_and_verify() {
    let d = Bn254::from_u32(D as u32);
    let air = EcAddAir { d };
    let trace = generate_trace_rows(d, &sample_pairs(8));

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
fn add_rejects_tampered_trace() {
    let d = Bn254::from_u32(D as u32);
    let air = EcAddAir { d };
    let mut trace = generate_trace_rows(d, &sample_pairs(8));

    // Corrupt the `x3` output of the first row.
    trace.values[4] += Bn254::ONE;

    let config = build_config(3);
    let _ = prove(&config, &air, trace, &[]);
}
