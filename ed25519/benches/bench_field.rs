//! Field-arithmetic microbenchmarks for the ed25519 base field (`p = 2^255 - 19`).
//!
//! Mirrors `bn254/benches/bench_field.rs` so the two 255-bit fields (and the
//! small 31/64-bit fields) are directly comparable under the same
//! `p3-field-testing` harness. Run with:
//!
//! ```text
//! CARGO_TARGET_DIR=/home/ubuntu/alp-target cargo bench -p p3-ed25519
//! ```
//!
//! `latency` = a dependent chain (one op waits on the previous — the worst case);
//! `throughput` = independent ops pipelined (what the prover actually sees).

use criterion::{Criterion, criterion_group, criterion_main};
use p3_ed25519::Ed25519;
use p3_field_testing::bench_func::{
    benchmark_add_latency, benchmark_add_throughput, benchmark_inv, benchmark_sqrt,
    benchmark_sub_latency, benchmark_sub_throughput,
};
use p3_field_testing::{benchmark_halve, benchmark_mul_latency, benchmark_mul_throughput};

type F = Ed25519;

fn bench_field(c: &mut Criterion) {
    let name = "Ed25519";
    const REPS: usize = 100;
    benchmark_halve::<F, REPS>(c, name);
    benchmark_inv::<F>(c, name);
    benchmark_sqrt::<F>(c, name);

    // Each throughput round runs 10 operations, so latency uses 10× the reps.
    const L_REPS: usize = 10 * REPS;
    benchmark_add_latency::<F, L_REPS>(c, name);
    benchmark_add_throughput::<F, REPS>(c, name);
    benchmark_sub_latency::<F, L_REPS>(c, name);
    benchmark_sub_throughput::<F, REPS>(c, name);
    benchmark_mul_latency::<F, L_REPS>(c, name);
    benchmark_mul_throughput::<F, REPS>(c, name);
}

criterion_group!(ed25519_arithmetic, bench_field);
criterion_main!(ed25519_arithmetic);
