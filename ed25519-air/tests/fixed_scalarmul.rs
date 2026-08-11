//! Structural tests for `EcFixedScalarMulAir` over `p3-bn254`'s field.
//!
//! See `tests/add.rs` for why `p3-bn254` and the `d` constant here are placeholders rather
//! than the real ed25519 base field / curve parameter.

use std::marker::PhantomData;

use p3_bn254::Bn254;
use p3_challenger::{CanObserve, CanSample, CanSampleBits, FieldChallenger};
use p3_commit::testing::TrivialPcs;
use p3_dft::Radix2DitParallel;
use p3_ed25519_air::fixed_scalarmul::{
    EcFixedScalarMulAir, NUM_EC_FIXED_SCALARMUL_COLS, generate_preprocessed_trace,
    generate_trace_rows,
};
use p3_field::{PrimeCharacteristicRing, PrimeField};
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

const NUM_ROWS: usize = 8;
const DEGREE_BITS: usize = 3; // log2(NUM_ROWS)

/// `scalar = 0b0110_1101 = 109`, least-significant bit first, 8 rows.
fn scalar_bits() -> Vec<bool> {
    let scalar = 0b0110_1101u8;
    (0..NUM_ROWS).map(|i| (scalar >> i) & 1 == 1).collect()
}

#[test]
fn fixed_scalarmul_prove_and_verify() {
    let d = Bn254::from_u32(D as u32);
    let base = (Bn254::from_u32(3), Bn254::from_u32(5));
    let air = EcFixedScalarMulAir {
        d,
        base,
        num_rows: NUM_ROWS,
    };
    let trace = generate_trace_rows(d, base, &scalar_bits());

    let config = build_config(DEGREE_BITS);
    let (preprocessed_prover_data, preprocessed_vk) =
        setup_preprocessed::<Config, _>(&config, &air, DEGREE_BITS).unwrap();

    let proof = prove_with_preprocessed(&config, &air, trace, &[], Some(&preprocessed_prover_data));
    verify_with_preprocessed(&config, &air, &proof, &[], Some(&preprocessed_vk))
        .expect("verification failed");
}

// `scalar` and the result are pure witness now (no public values); tampering with either is
// just another way of tampering with the trace, caught the same way
// `fixed_scalarmul_rejects_tampered_bit` below catches a corrupted `bit` column. In debug
// builds, `prove_with_preprocessed` itself checks constraints and panics on a tampered trace,
// before a bad proof could even reach `verify_with_preprocessed`; see the comment on
// `add_rejects_tampered_trace` in `tests/add.rs`.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "constraints not satisfied on row")]
fn fixed_scalarmul_rejects_tampered_result() {
    let d = Bn254::from_u32(D as u32);
    let base = (Bn254::from_u32(3), Bn254::from_u32(5));
    let air = EcFixedScalarMulAir {
        d,
        base,
        num_rows: NUM_ROWS,
    };
    let mut trace = generate_trace_rows(d, base, &scalar_bits());

    // Corrupt the last row's `result_x` column (index 10 within the row), desyncing it from
    // the mux expression it's constrained equal to.
    let last_row_result_x = (NUM_ROWS - 1) * NUM_EC_FIXED_SCALARMUL_COLS + 10;
    trace.values[last_row_result_x] += Bn254::ONE;

    let config = build_config(DEGREE_BITS);
    let (preprocessed_prover_data, _) =
        setup_preprocessed::<Config, _>(&config, &air, DEGREE_BITS).unwrap();
    let _ = prove_with_preprocessed(&config, &air, trace, &[], Some(&preprocessed_prover_data));
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "constraints not satisfied on row")]
fn fixed_scalarmul_rejects_tampered_bit() {
    let d = Bn254::from_u32(D as u32);
    let base = (Bn254::from_u32(3), Bn254::from_u32(5));
    let air = EcFixedScalarMulAir {
        d,
        base,
        num_rows: NUM_ROWS,
    };
    let mut trace = generate_trace_rows(d, base, &scalar_bits());

    // Flip the `bit` column of the first row; this should desync the boolean/mux constraints
    // (or the scalar-reconstruction sum) from the rest of the trace.
    trace.values[0] = Bn254::ONE - trace.values[0];

    let config = build_config(DEGREE_BITS);
    let (preprocessed_prover_data, _) =
        setup_preprocessed::<Config, _>(&config, &air, DEGREE_BITS).unwrap();
    let _ = prove_with_preprocessed(&config, &air, trace, &[], Some(&preprocessed_prover_data));
}

/// Cheap regression guard against the claimed witness-width reduction (12 columns, down from
/// `EcScalarMulAir`'s 23) silently regressing.
#[test]
fn fixed_scalarmul_witness_width_is_12() {
    assert_eq!(NUM_EC_FIXED_SCALARMUL_COLS, 12);
}

#[test]
fn fixed_scalarmul_preprocessed_trace_matches_doubling_ladder() {
    let d = Bn254::from_u32(D as u32);
    let base = (Bn254::from_u32(3), Bn254::from_u32(5));

    let preprocessed = generate_preprocessed_trace(d, base, NUM_ROWS);

    let mut expected = Vec::with_capacity(NUM_ROWS);
    let (mut x, mut y) = base;
    for _ in 0..NUM_ROWS {
        expected.push((x, y));
        let (nx, ny) = p3_ed25519_air::add::add_points(d, (x, y), (x, y));
        x = nx;
        y = ny;
    }

    for (row, &(ex, ey)) in expected.iter().enumerate() {
        assert_eq!(
            preprocessed.values[row * 2],
            ex,
            "x_dbl mismatch at row {row}"
        );
        assert_eq!(
            preprocessed.values[row * 2 + 1],
            ey,
            "y_dbl mismatch at row {row}"
        );
    }
}
