//! QA encoding throughput over ed25519, and the same message length through
//! the two codes already in the tree (Brakedown, Reed–Solomon).
//!
//! The claim under test is the paper's Table 1: QA encoding at rate 1/2 costs
//! `N` multiplications and `2N log N` additions, which should make it the
//! cheapest of the three per encoded element over an expensive 255-bit field.
//!
//! Run: `cargo bench -p p3-qa-code --bench qa_encode --features parallel`

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use p3_brakedown::{BrakedownCode, LinearCode, ReedSolomonCode};
use p3_ed25519::Ed25519;
use p3_matrix::dense::RowMajorMatrix;
use p3_qa_code::QaCode;
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;

type F = Ed25519;
type EF = Ed25519;

/// Matrix rows encoded together (the tensor PCS width).
const WIDTH: usize = 8;

fn rand_matrix(rng: &mut ChaCha8Rng, height: usize, width: usize) -> RowMajorMatrix<F> {
    RowMajorMatrix::new((0..height * width).map(|_| rng.random::<F>()).collect(), width)
}

fn bench_qa(c: &mut Criterion) {
    let mut group = c.benchmark_group("qa_encode_ed25519");
    let mut rng = ChaCha8Rng::seed_from_u64(0xE1);

    for log_n in [10, 12, 14, 16] {
        let n = 1usize << log_n;
        group.throughput(Throughput::Elements((n * WIDTH) as u64));

        for index in [2, 4] {
            let code = QaCode::<F>::new_random(log_n, index, &mut ChaCha8Rng::seed_from_u64(1));
            let msg = rand_matrix(&mut rng, n, WIDTH);
            group.bench_with_input(
                BenchmarkId::new(format!("qa_rate_1_over_{index}"), log_n),
                &log_n,
                |b, _| b.iter(|| code.encode(&msg)),
            );
        }
    }
    group.finish();
}

/// The other two codes at the same message length and rate 1/2-ish, so the
/// per-element cost is directly comparable.
fn bench_alternatives(c: &mut Criterion) {
    let mut group = c.benchmark_group("code_encode_ed25519_comparison");
    let mut rng = ChaCha8Rng::seed_from_u64(0xE2);

    for log_n in [10, 12] {
        let n = 1usize << log_n;
        let msg = rand_matrix(&mut rng, n, WIDTH);
        group.throughput(Throughput::Elements((n * WIDTH) as u64));

        let qa = QaCode::<F>::new_random(log_n, 2, &mut ChaCha8Rng::seed_from_u64(1));
        group.bench_with_input(BenchmarkId::new("qa_rate_1_over_2", log_n), &log_n, |b, _| {
            b.iter(|| qa.encode(&msg))
        });

        let brakedown =
            BrakedownCode::<F, EF>::from_message_len(n, &mut ChaCha8Rng::seed_from_u64(2));
        group.bench_with_input(BenchmarkId::new("brakedown", log_n), &log_n, |b, _| {
            b.iter(|| brakedown.encode(&msg))
        });

        // Reed–Solomon rate 1/2 — superlinear (Vandermonde), included as the
        // reference the tensor PCS already ships (Ligero).
        let rs = ReedSolomonCode::<F>::new(n, 2 * n);
        group.bench_with_input(BenchmarkId::new("reed_solomon_rate_1_over_2", log_n), &log_n, |b, _| {
            b.iter(|| rs.encode(&msg))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_qa, bench_alternatives);
criterion_main!(benches);
