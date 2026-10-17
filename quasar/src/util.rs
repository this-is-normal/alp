//! Multilinear helpers for the code-switching layer.
//!
//! # Index convention
//!
//! Everything here is **MSB-first**: for an `n`-variable table, variable `k`
//! sits at bit `n−1−k` of the array index, so variable `0` is the most
//! significant. This matches `Poly::fix_prefix_var_mut` (which binds variable
//! `0` first), `p3-basefold`'s `bind_msb`, and the row/column split of the
//! committed matrix (`M[i][j] = evals[i·N₂ + j]`, so the `n₁` row variables are
//! the leading ones).

use alloc::vec::Vec;

use p3_field::Field;
use p3_maybe_rayon::prelude::*;

/// `eq(a, x) = Π_k (a_k·x_k + (1−a_k)(1−x_k))` over the hypercube, MSB-first.
pub fn eq_table<F: Field>(a: &[F]) -> Vec<F> {
    tensor_table(a, |ak| (F::ONE - ak, ak))
}

/// `eq(a, b)` at two arbitrary points.
pub fn eq_eval<F: Field>(a: &[F], b: &[F]) -> F {
    debug_assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(&x, &y)| x * y + (F::ONE - x) * (F::ONE - y))
        .product()
}

/// The multilinear extension of the Walsh–Hadamard matrix,
/// `H̃_n(a, b) = Π_k (1 − 2·a_k·b_k)`.
///
/// `H_n = H^{⊗n}` has entries `H_n[i][j] = (−1)^{popcount(i & j)}`, whose
/// unique multilinear extension factorizes because `H`'s extension is
/// `1 − 2ab`. Costs `O(n)` — this is what lets the Quasar verifier check a WHT
/// relation in logarithmic time. Note `popcount(i & j)` is invariant under
/// permuting bit positions, so `H̃` needs no index convention of its own.
pub fn wht_mle_eval<F: Field>(a: &[F], b: &[F]) -> F {
    debug_assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(&x, &y)| F::ONE - F::TWO * x * y)
        .product()
}

/// The table `[H̃_n(a, x)]` for all `x` in the hypercube, MSB-first.
///
/// At a boolean `x` the factor is `1` when `x_k = 0` and `1 − 2a_k` when
/// `x_k = 1`.
pub fn wht_mle_table<F: Field>(a: &[F]) -> Vec<F> {
    tensor_table(a, |ak| (F::ONE, F::ONE - F::TWO * ak))
}

/// Build `Π_k factor(a_k)[x_k]` over the hypercube, MSB-first.
///
/// Variables are consumed in order, each becoming the new least significant
/// bit, so the first variable ends up most significant.
fn tensor_table<F: Field, G: Fn(F) -> (F, F)>(a: &[F], factor: G) -> Vec<F> {
    let mut table = alloc::vec![F::ONE];
    for &ak in a {
        let (f0, f1) = factor(ak);
        // Each existing entry splits into an adjacent pair, independently.
        let mut next = F::zero_vec(2 * table.len());
        next.par_chunks_exact_mut(2)
            .zip(table.par_iter())
            .for_each(|(pair, &t)| {
                pair[0] = t * f0;
                pair[1] = t * f1;
            });
        table = next;
    }
    table
}

/// `f̃(z)` for a multilinear given by its hypercube evaluations.
pub fn eval_mle<F: Field>(evals: &[F], z: &[F]) -> F {
    debug_assert_eq!(evals.len(), 1 << z.len());
    eq_table(z)
        .par_iter()
        .zip(evals.par_iter())
        .map(|(&e, &f)| e * f)
        .sum()
}

/// `Σ_i w_i · m_i`, the standard inner product.
pub fn dot<F: Field>(a: &[F], b: &[F]) -> F {
    debug_assert_eq!(a.len(), b.len());
    // Left serial: the Quasar verifier calls this once per opened column with
    // `N₁` (tens to hundreds) elements, where a task would cost more than the
    // work, and the *caller* parallelises across columns instead.
    a.iter().zip(b).map(|(&x, &y)| x * y).sum()
}

/// `[1, x, x², …, x^{n−1}]`.
pub fn powers<F: Field>(x: F, n: usize) -> Vec<F> {
    let mut out = Vec::with_capacity(n);
    let mut p = F::ONE;
    for _ in 0..n {
        out.push(p);
        p *= x;
    }
    out
}

/// `eq(⟨index⟩, point)` for a hypercube point given as an integer, MSB-first.
///
/// Used to evaluate the sparse selector polynomial without materializing it.
pub fn eq_at_index<F: Field>(index: usize, num_vars: usize, point: &[F]) -> F {
    debug_assert_eq!(point.len(), num_vars);
    point
        .iter()
        .enumerate()
        .map(|(k, &pk)| {
            if (index >> (num_vars - 1 - k)) & 1 == 1 {
                pk
            } else {
                F::ONE - pk
            }
        })
        .product()
}

#[cfg(test)]
mod tests {
    use p3_ed25519::Ed25519;
    use p3_field::PrimeCharacteristicRing;
    use p3_multilinear_util::point::Point;
    use p3_multilinear_util::poly::Poly;
    use p3_qa_code::wht;
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    use super::*;

    type F = Ed25519;

    fn rand_vec(rng: &mut ChaCha8Rng, n: usize) -> Vec<F> {
        (0..n).map(|_| rng.random::<F>()).collect()
    }

    /// Our MSB-first `eq_table` must agree with the shared `Poly`/`Point`
    /// helpers, or the row/column split and the Basefold openings would use
    /// different variable orders.
    #[test]
    fn eq_table_matches_shared_helpers() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for n in 1..8 {
            let z = rand_vec(&mut rng, n);
            assert_eq!(
                eq_table(&z),
                Poly::<F>::new_from_point(&z, F::ONE).as_slice().to_vec(),
                "n={n}"
            );

            let w = rand_vec(&mut rng, n);
            assert_eq!(eq_eval(&z, &w), Point::<F>::eval_eq(&z, &w), "n={n}");
            assert_eq!(eval_mle(&eq_table(&z), &w), eq_eval(&z, &w), "n={n}");
        }
    }

    /// `H̃` must be the multilinear extension of the WHT matrix: agree with
    /// `(−1)^{popcount(i & j)}` at boolean points, and with the table form.
    #[test]
    fn wht_mle_extends_the_wht_matrix() {
        let mut rng = ChaCha8Rng::seed_from_u64(1);
        for n in 1..7 {
            let len = 1usize << n;
            // Boolean agreement.
            for i in 0..len {
                let a: Vec<F> = (0..n)
                    .map(|k| F::from_bool((i >> (n - 1 - k)) & 1 == 1))
                    .collect();
                for j in 0..len {
                    let b: Vec<F> = (0..n)
                        .map(|k| F::from_bool((j >> (n - 1 - k)) & 1 == 1))
                        .collect();
                    let expected = if (i & j).count_ones() % 2 == 0 {
                        F::ONE
                    } else {
                        F::NEG_ONE
                    };
                    assert_eq!(wht_mle_eval(&a, &b), expected, "n={n} i={i} j={j}");
                }
            }
            // Table form and arbitrary-point form agree.
            let a = rand_vec(&mut rng, n);
            let b = rand_vec(&mut rng, n);
            assert_eq!(eval_mle(&wht_mle_table(&a), &b), wht_mle_eval(&a, &b));
        }
    }

    /// The identity the encoding sumcheck rests on:
    /// `WHT(v)~(γ) = Σ_x H̃(γ, x)·ṽ(x)`.
    #[test]
    fn wht_relation_holds_as_a_sumcheck_claim() {
        let mut rng = ChaCha8Rng::seed_from_u64(2);
        for n in 1..8 {
            let v = rand_vec(&mut rng, 1 << n);
            let mut transformed = v.clone();
            wht(&mut transformed);

            let gamma = rand_vec(&mut rng, n);
            let lhs = eval_mle(&transformed, &gamma);
            let rhs = dot(&wht_mle_table(&gamma), &v);
            assert_eq!(lhs, rhs, "n={n}");
        }
    }

    #[test]
    fn eq_at_index_matches_the_table() {
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        let n = 6;
        let z = rand_vec(&mut rng, n);
        let table = eq_table(&z);
        for (i, &t) in table.iter().enumerate() {
            assert_eq!(eq_at_index(i, n, &z), t, "index {i}");
        }
    }
}
