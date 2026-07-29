//! Structural tests for `EcScalarMulAir` over `p3-bn254`'s field.
//!
//! See `tests/add.rs` for why `p3-bn254` and the `d` constant here are placeholders rather
//! than the real ed25519 base field / curve parameter.

use std::marker::PhantomData;

use p3_bn254::Bn254;
use p3_challenger::{CanObserve, CanSample, CanSampleBits, FieldChallenger};
use p3_commit::testing::TrivialPcs;
use p3_dft::Radix2DitParallel;
use p3_ed25519_air::scalarmul::{EcScalarMulAir, generate_trace_rows};
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

/// `scalar = 0b0110_1101 = 109`, least-significant bit first, 8 rows.
fn scalar_bits() -> Vec<bool> {
    let scalar = 0b0110_1101u8;
    (0..8).map(|i| (scalar >> i) & 1 == 1).collect()
}

#[test]
fn scalarmul_prove_and_verify() {
    let d = Bn254::from_u32(D as u32);
    let air = EcScalarMulAir { d };
    let base = (Bn254::from_u32(3), Bn254::from_u32(5));
    let result = generate_trace_rows(d, base, &scalar_bits());

    let config = build_config(3);
    let proof = prove(&config, &air, result.trace, &result.public_values);
    verify(&config, &air, &proof, &result.public_values).expect("verification failed");
}

// `TrivialPcs` (a toy PCS for structural tests, not a real commitment scheme) checks
// opening consistency with a plain `assert_eq!` rather than returning a `Result::Err`, so a
// mismatched public value surfaces here as a panic rather than an `Err` from `verify`.
#[test]
#[should_panic]
fn scalarmul_rejects_wrong_public_result() {
    let d = Bn254::from_u32(D as u32);
    let air = EcScalarMulAir { d };
    let base = (Bn254::from_u32(3), Bn254::from_u32(5));
    let result = generate_trace_rows(d, base, &scalar_bits());

    let config = build_config(3);
    let proof = prove(&config, &air, result.trace, &result.public_values);

    let mut wrong_public_values = result.public_values;
    wrong_public_values[3] += Bn254::ONE; // corrupt the claimed result_x
    let _ = verify(&config, &air, &proof, &wrong_public_values);
}

// In debug builds, `prove` itself checks constraints and panics on a tampered trace,
// before a bad proof could even reach `verify`; see the comment on `add_rejects_tampered_trace`
// in `tests/add.rs`.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "constraints not satisfied on row")]
fn scalarmul_rejects_tampered_bit() {
    let d = Bn254::from_u32(D as u32);
    let air = EcScalarMulAir { d };
    let base = (Bn254::from_u32(3), Bn254::from_u32(5));
    let result = generate_trace_rows(d, base, &scalar_bits());
    let mut trace = result.trace;

    // Flip the `bit` column of the first row; this should desync the boolean/mux constraints
    // (or the scalar-reconstruction sum) from the public values.
    trace.values[0] = Bn254::ONE - trace.values[0];

    let config = build_config(3);
    let _ = prove(&config, &air, trace, &result.public_values);
}
