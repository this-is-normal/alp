//! The Walsh–Hadamard transform, batched over the columns of a row-major
//! matrix.
//!
//! `WHT(v) = H^{⊗n} · v` with `H = [[1, 1], [1, -1]]`. The Sylvester–Hadamard
//! matrix is symmetric and the decimation-in-time butterfly below leaves the
//! output in natural order — there is no bit-reversal permutation, unlike an
//! FFT. Additions only: no field multiplications and no twiddle factors, which
//! is why QA encoding is cheap.
//!
//! Both the vector and the batched form are unnormalized: `WHT(WHT(v)) = 2ⁿ·v`.
//! Callers avoid the inverse transform by folding `2^{-n}` into the public QA
//! coefficient vectors (see [`crate::QaCode`]).

use p3_field::Algebra;
use p3_field::Field;
use p3_maybe_rayon::prelude::*;

/// Blocks per level below which we parallelize *inside* a block instead of
/// across blocks. The last few levels have very few blocks but very long
/// halves, so the two strategies cover opposite ends of the transform.
const MIN_BLOCKS_FOR_OUTER_PAR: usize = 32;

/// Half-length below which a block's butterfly runs serially: rayon overhead
/// dominates for short halves.
const MIN_HALF_FOR_INNER_PAR: usize = 1 << 10;

/// In-place WHT of a single vector. `x.len()` must be a power of two.
pub fn wht<F: Field>(x: &mut [F]) {
    wht_batched(x, 1);
}

/// In-place WHT of every column of a row-major `height × width` matrix stored
/// in `values` (so `values.len() == height * width`), where `height` is a power
/// of two.
///
/// This is the form the code needs: the [`LinearCode`](p3_brakedown::LinearCode)
/// convention makes each *column* an independent message, so the butterflies
/// pair rows and stride by `width`.
pub fn wht_batched<A>(values: &mut [A], width: usize)
where
    A: Copy + Send + Sync + core::ops::Add<Output = A> + core::ops::Sub<Output = A>,
{
    assert!(width > 0, "width must be positive");
    assert_eq!(values.len() % width, 0, "values must be a whole matrix");
    let height = values.len() / width;
    assert!(
        height.is_power_of_two(),
        "WHT height must be a power of two, got {height}"
    );

    // `half` = elements in one half of a butterfly block = (rows per half)·width.
    let mut half = width;
    while 2 * half <= values.len() {
        let block = 2 * half;
        let blocks = values.len() / block;

        if blocks >= MIN_BLOCKS_FOR_OUTER_PAR {
            values
                .par_chunks_exact_mut(block)
                .for_each(|blk| butterfly_serial(blk, half));
        } else if half >= MIN_HALF_FOR_INNER_PAR {
            values
                .chunks_exact_mut(block)
                .for_each(|blk| butterfly_parallel(blk, half));
        } else {
            values
                .chunks_exact_mut(block)
                .for_each(|blk| butterfly_serial(blk, half));
        }

        half = block;
    }
}

/// `(lo, hi) ← (lo + hi, lo − hi)` over one butterfly block.
#[inline]
fn butterfly_serial<A>(block: &mut [A], half: usize)
where
    A: Copy + core::ops::Add<Output = A> + core::ops::Sub<Output = A>,
{
    let (lo, hi) = block.split_at_mut(half);
    for (l, r) in lo.iter_mut().zip(hi.iter_mut()) {
        let (u, v) = (*l, *r);
        *l = u + v;
        *r = u - v;
    }
}

#[inline]
fn butterfly_parallel<A>(block: &mut [A], half: usize)
where
    A: Copy + Send + Sync + core::ops::Add<Output = A> + core::ops::Sub<Output = A>,
{
    let (lo, hi) = block.split_at_mut(half);
    lo.par_iter_mut().zip(hi.par_iter_mut()).for_each(|(l, r)| {
        let (u, v) = (*l, *r);
        *l = u + v;
        *r = u - v;
    });
}

/// The naive `H^{⊗n}·v`, for tests: `out[i] = Σ_j (−1)^{popcount(i & j)} v[j]`.
///
/// Quadratic; only for small instances.
pub fn wht_naive<F: Field, A: Algebra<F> + Copy>(v: &[A]) -> alloc::vec::Vec<A> {
    let n = v.len();
    assert!(n.is_power_of_two());
    (0..n)
        .map(|i| {
            v.iter().enumerate().fold(A::ZERO, |acc, (j, &x)| {
                if (i & j).count_ones() % 2 == 0 {
                    acc + x
                } else {
                    acc - x
                }
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use p3_baby_bear::BabyBear;
    use p3_field::PrimeCharacteristicRing;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    use super::*;

    type F = BabyBear;

    fn random_vec(rng: &mut ChaCha8Rng, n: usize) -> Vec<F> {
        use rand::RngExt;
        (0..n).map(|_| rng.random()).collect()
    }

    #[test]
    fn matches_naive() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for log_n in 0..8 {
            let v = random_vec(&mut rng, 1 << log_n);
            let expected = wht_naive::<F, F>(&v);
            let mut got = v.clone();
            wht(&mut got);
            assert_eq!(got, expected, "log_n = {log_n}");
        }
    }

    #[test]
    fn involution_up_to_scale() {
        let mut rng = ChaCha8Rng::seed_from_u64(1);
        for log_n in 0..10 {
            let n = 1 << log_n;
            let v = random_vec(&mut rng, n);
            let mut got = v.clone();
            wht(&mut got);
            wht(&mut got);
            let scale = F::from_usize(n);
            for (g, o) in got.iter().zip(&v) {
                assert_eq!(*g, *o * scale);
            }
        }
    }

    /// The batched form must transform each column independently.
    #[test]
    fn batched_matches_per_column() {
        let mut rng = ChaCha8Rng::seed_from_u64(2);
        for (log_h, width) in [(1, 1), (3, 5), (6, 3), (10, 7)] {
            let height = 1usize << log_h;
            let values = random_vec(&mut rng, height * width);

            let mut got = values.clone();
            wht_batched(&mut got, width);

            for c in 0..width {
                let column: Vec<F> = (0..height).map(|r| values[r * width + c]).collect();
                let mut expected = column;
                wht(&mut expected);
                for (r, e) in expected.iter().enumerate() {
                    assert_eq!(got[r * width + c], *e, "h={height} w={width} col={c}");
                }
            }
        }
    }

    /// `WHT` is `F`-linear, which is what makes the QA code a linear code.
    #[test]
    fn is_linear() {
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        let n = 64;
        let (a, b) = (random_vec(&mut rng, n), random_vec(&mut rng, n));
        let (mut ta, mut tb) = (a.clone(), b.clone());
        wht(&mut ta);
        wht(&mut tb);

        let mut sum: Vec<F> = a.iter().zip(&b).map(|(x, y)| *x + *y).collect();
        wht(&mut sum);
        for (s, (x, y)) in sum.iter().zip(ta.iter().zip(&tb)) {
            assert_eq!(*s, *x + *y);
        }
    }
}
