//! Per-primitive ed25519 operation counts, so the cost model's arithmetic terms
//! are measured constants rather than fitted ones.
//!
//! Usage: `cargo run -p pcs-profile --release --bin ed25519_primitives`
//! (needs the default `count-ops` feature).

use p3_ed25519::Ed25519;
use p3_field::{Field, PrimeCharacteristicRing, batch_multiplicative_inverse};
use p3_baby_bear::{BabyBear, default_babybear_poseidon2_16};
use p3_stark_stats::{CountingPerm, HashCounts};
use p3_symmetric::{CryptographicHasher, PaddingFreeSponge, SplitFieldHasher};
use pcs_profile::{BB_DIGEST, BB_RATE, BB_WIDTH, reset_arith, snapshot_arith};
use p3_quasar::{eq_eval, eq_table, wht_mle_eval};
use p3_sumcheck::generic_degree::RoundPolyInterpolator;
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;

fn measure<T>(label: &str, reps: usize, mut f: impl FnMut() -> T) {
    reset_arith();
    for _ in 0..reps {
        core::hint::black_box(f());
    }
    let (mul, sq, inv, add) = snapshot_arith();
    let r = reps as f64;
    println!(
        "{label:<44} mul {:>9.2}  sq {:>8.2}  inv {:>6.2}  add {:>9.2}   (per call, {reps} reps)",
        mul as f64 / r,
        sq as f64 / r,
        inv as f64 / r,
        add as f64 / r
    );
}

fn main() {
    let mut rng = ChaCha8Rng::seed_from_u64(7);
    let a: Ed25519 = rng.random();
    let b: Ed25519 = rng.random();
    let vec3: Vec<Ed25519> = (0..3).map(|_| rng.random()).collect();
    let vec4: Vec<Ed25519> = (0..4).map(|_| rng.random()).collect();

    println!("# ed25519 primitive costs (counters ON)");
    measure("mul", 1000, || a * b);
    measure("square", 1000, || a.square());
    measure("inverse (Fermat a^(p-2))", 100, || a.inverse());
    measure("batch_multiplicative_inverse(3)", 100, || {
        batch_multiplicative_inverse(&vec3)
    });
    measure("batch_multiplicative_inverse(4)", 100, || {
        batch_multiplicative_inverse(&vec4)
    });

    // Sumcheck-round interpolation: the verifier's dominant non-column cost.
    for degree in [2usize, 3] {
        let label = format!("RoundPolyInterpolator::new({degree})");
        measure(&label, 100, || RoundPolyInterpolator::<Ed25519>::new(degree));

        let interp = RoundPolyInterpolator::<Ed25519>::new(degree);
        let evals: Vec<Ed25519> = (0..degree).map(|_| rng.random()).collect();
        let label = format!("RoundPolyInterpolator::eval(degree {degree})");
        measure(&label, 100, || interp.eval(&evals, a, b));
    }

    // Multilinear helpers the Quasar verifier calls.
    for n in [8usize, 16] {
        let z: Vec<Ed25519> = (0..n).map(|_| rng.random()).collect();
        let w: Vec<Ed25519> = (0..n).map(|_| rng.random()).collect();
        measure(&format!("eq_table({n} vars)"), 50, || eq_table(&z));
        measure(&format!("eq_eval({n} vars)"), 200, || eq_eval(&z, &w));
        measure(&format!("wht_mle_eval({n} vars)"), 200, || wht_mle_eval(&z, &w));
    }

    // Leaf hashing: does the Montgomery -> canonical limb split cost field muls?
    let counts = HashCounts::new();
    let hasher: SplitFieldHasher<
        Ed25519,
        BabyBear,
        PaddingFreeSponge<CountingPerm<_>, BB_WIDTH, BB_RATE, BB_DIGEST>,
    > = SplitFieldHasher::new(PaddingFreeSponge::new(CountingPerm::with_counter(
        default_babybear_poseidon2_16(),
        counts.leaf_hash.clone(),
    )));
    for m in [1usize, 16, 64, 256] {
        let row: Vec<Ed25519> = (0..m).map(|_| rng.random()).collect();
        let label = format!("leaf hash of {m} elements");
        measure(&label, 20, || {
            let out: [BabyBear; BB_DIGEST] = hasher.hash_iter(row.iter().copied());
            out
        });
    }
}
