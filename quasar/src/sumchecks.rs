//! The two sumcheck shapes Quasar needs, as [`RoundProver`] states driven by
//! the shared generic-degree engine in `p3-sumcheck`.
//!
//! - [`Deg2Prover`]: `Σ_x Σ_p a_p(x)·b_p(x)` — the selector check, the inner
//!   WHT layer, and the claim-reduction batcher.
//! - [`Deg3Prover`]: `Σ_x Σ_p a_p(x)·b_p(x)·c_p(x)` — the outer WHT layer,
//!   whose three factors are `H̃(γ,·)`, the batched public coefficients, and the
//!   encoding witness `WHT(v)`.
//!
//! Both bind variable `0` first (MSB-first, matching [`crate::util`]), and
//! report `h(0), h(2), …, h(degree)` as the engine expects — `h(1)` is derived
//! by the verifier from the running claim.
//!
//! Each round is a single parallel reduction over the two halves, computing
//! every reported node at once: the naive shape (interpolate each factor at
//! each node into a fresh `Vec`, then zip) reads the same `O(2^k)` data once
//! per node and allocates per factor per node.

use alloc::vec::Vec;

use p3_field::Field;
use p3_maybe_rayon::prelude::*;
use p3_multilinear_util::poly::Poly;
use p3_sumcheck::generic_degree::RoundProver;

/// Additive accumulator for a multi-value reduction.
///
/// `p3-maybe-rayon`'s two backends disagree on `reduce`: rayon's takes an
/// identity plus a combiner, `core::iter`'s takes only a combiner. `sum` has the
/// same shape in both, so accumulating through this type is what lets one
/// expression compile serial and parallel.
#[derive(Copy, Clone)]
struct Sums<F, const N: usize>([F; N]);

impl<F: Field, const N: usize> Default for Sums<F, N> {
    fn default() -> Self {
        Self([F::ZERO; N])
    }
}

impl<F: Field, const N: usize> core::ops::Add for Sums<F, N> {
    type Output = Self;

    fn add(mut self, other: Self) -> Self {
        for (a, b) in self.0.iter_mut().zip(other.0) {
            *a += b;
        }
        self
    }
}

impl<F: Field, const N: usize> core::iter::Sum for Sums<F, N> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), |a, b| a + b)
    }
}

/// Multilinear interpolation of a factor's two halves at an integer node:
/// `f(t) = f(0) + (f(1) − f(0))·t`, so `t = 0` is the low half, `t = 2` is
/// `2f₁ − f₀`, and `t = 3` is `3f₁ − 2f₀`.
#[inline]
fn at_node<F: Field>(f0: F, f1: F, node: Node) -> F {
    match node {
        Node::Two => F::TWO * f1 - f0,
        Node::Three => F::from_u8(3) * f1 - F::TWO * f0,
    }
}

/// The off-`{0,1}` nodes a round polynomial is reported at. `h(0)` is the low
/// half verbatim and needs no interpolation; `h(1)` is never sent — the verifier
/// derives it from the running claim.
#[derive(Copy, Clone)]
enum Node {
    Two,
    Three,
}

/// `Σ_x Σ_p a_p(x)·b_p(x)` over the hypercube.
pub struct Deg2Prover<F: Field> {
    pairs: Vec<(Poly<F>, Poly<F>)>,
}

impl<F: Field> Deg2Prover<F> {
    /// Every pair must have the same number of variables.
    pub fn new(pairs: Vec<(Poly<F>, Poly<F>)>) -> Self {
        assert!(!pairs.is_empty(), "sumcheck needs at least one term");
        let n = pairs[0].0.num_evals();
        assert!(
            pairs
                .iter()
                .all(|(a, b)| a.num_evals() == n && b.num_evals() == n),
            "all factors must share the variable space"
        );
        Self { pairs }
    }

    /// The claimed sum, for the caller to put in the transcript.
    pub fn claimed_sum(&self) -> F {
        self.pairs
            .iter()
            .map(|(a, b)| {
                a.as_slice()
                    .par_iter()
                    .zip(b.as_slice().par_iter())
                    .map(|(&x, &y)| x * y)
                    .sum::<F>()
            })
            .sum()
    }

    /// Both reported evaluations in ONE pass over the halves. Evaluating each
    /// node separately would read the same `O(2^k)` data twice and allocate an
    /// interpolated copy of every factor per node.
    fn round_evals(&self) -> [F; 2] {
        self.pairs
            .iter()
            .map(|(a, b)| {
                let half = a.num_evals() / 2;
                let (a0, a1) = a.as_slice().split_at(half);
                let (b0, b1) = b.as_slice().split_at(half);
                a0.par_iter()
                    .zip(a1.par_iter())
                    .zip(b0.par_iter().zip(b1.par_iter()))
                    .map(|((&x0, &x1), (&y0, &y1))| {
                        Sums([
                            x0 * y0,
                            at_node(x0, x1, Node::Two) * at_node(y0, y1, Node::Two),
                        ])
                    })
                    .sum::<Sums<F, 2>>()
            })
            .fold(Sums::default(), |a: Sums<F, 2>, b| a + b)
            .0
    }
}

impl<F: Field> RoundProver<F> for Deg2Prover<F> {
    fn fold(&mut self, r: F) {
        for (a, b) in &mut self.pairs {
            a.fix_prefix_var_mut(r);
            b.fix_prefix_var_mut(r);
        }
    }

    fn round_poly(&self) -> Vec<F> {
        self.round_evals().to_vec()
    }
}

/// `Σ_x Σ_p a_p(x)·b_p(x)·c_p(x)` over the hypercube.
pub struct Deg3Prover<F: Field> {
    triples: Vec<(Poly<F>, Poly<F>, Poly<F>)>,
}

impl<F: Field> Deg3Prover<F> {
    pub fn new(triples: Vec<(Poly<F>, Poly<F>, Poly<F>)>) -> Self {
        assert!(!triples.is_empty(), "sumcheck needs at least one term");
        let n = triples[0].0.num_evals();
        assert!(
            triples
                .iter()
                .all(|(a, b, c)| a.num_evals() == n && b.num_evals() == n && c.num_evals() == n),
            "all factors must share the variable space"
        );
        Self { triples }
    }

    pub fn claimed_sum(&self) -> F {
        self.triples
            .iter()
            .map(|(a, b, c)| {
                a.as_slice()
                    .par_iter()
                    .zip(b.as_slice().par_iter())
                    .zip(c.as_slice().par_iter())
                    .map(|((&x, &y), &z)| x * y * z)
                    .sum::<F>()
            })
            .sum()
    }

    /// All three reported evaluations in one pass — see [`Deg2Prover::round_evals`].
    fn round_evals(&self) -> [F; 3] {
        self.triples
            .iter()
            .map(|(a, b, c)| {
                let half = a.num_evals() / 2;
                let (a0, a1) = a.as_slice().split_at(half);
                let (b0, b1) = b.as_slice().split_at(half);
                let (c0, c1) = c.as_slice().split_at(half);
                a0.par_iter()
                    .zip(a1.par_iter())
                    .zip(b0.par_iter().zip(b1.par_iter()))
                    .zip(c0.par_iter().zip(c1.par_iter()))
                    .map(|(((&x0, &x1), (&y0, &y1)), (&z0, &z1))| {
                        Sums([
                            x0 * y0 * z0,
                            at_node(x0, x1, Node::Two)
                                * at_node(y0, y1, Node::Two)
                                * at_node(z0, z1, Node::Two),
                            at_node(x0, x1, Node::Three)
                                * at_node(y0, y1, Node::Three)
                                * at_node(z0, z1, Node::Three),
                        ])
                    })
                    .sum::<Sums<F, 3>>()
            })
            .fold(Sums::default(), |a: Sums<F, 3>, b| a + b)
            .0
    }
}

impl<F: Field> RoundProver<F> for Deg3Prover<F> {
    fn fold(&mut self, r: F) {
        for (a, b, c) in &mut self.triples {
            a.fix_prefix_var_mut(r);
            b.fix_prefix_var_mut(r);
            c.fix_prefix_var_mut(r);
        }
    }

    fn round_poly(&self) -> Vec<F> {
        self.round_evals().to_vec()
    }
}

#[cfg(test)]
mod tests {
    use p3_challenger::{DuplexChallengerPF, FieldChallenger};
    use p3_ed25519::{
        ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519,
        Poseidon2Ed25519,
    };
    use p3_field::PrimeCharacteristicRing;
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    use super::*;
    use crate::util::{eval_mle, powers};

    type F = Ed25519;
    const WIDTH: usize = 3;
    type Perm = Poseidon2Ed25519<WIDTH>;
    type Ch = DuplexChallengerPF<Ed25519, Perm, WIDTH, 2>;

    fn challenger() -> Ch {
        let mut rng = ChaCha8Rng::seed_from_u64(0xC0FFEE);
        Ch::new(Poseidon2Ed25519::<WIDTH>::new_from_rng(
            ED25519_POSEIDON2_FULL_ROUNDS,
            ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
            &mut rng,
        ))
    }

    fn rand_poly(rng: &mut ChaCha8Rng, n: usize) -> Poly<F> {
        Poly::new((0..(1 << n)).map(|_| rng.random::<F>()).collect())
    }

    /// A degree-2 sumcheck must reduce to the product of the factors at the
    /// challenge point.
    #[test]
    fn deg2_roundtrip() {
        let mut rng = ChaCha8Rng::seed_from_u64(1);
        let n = 6;
        let pairs: Vec<(Poly<F>, Poly<F>)> = (0..3)
            .map(|_| (rand_poly(&mut rng, n), rand_poly(&mut rng, n)))
            .collect();
        let originals = pairs.clone();

        let mut prover = Deg2Prover::new(pairs);
        let claimed = prover.claimed_sum();

        let mut ch = challenger();
        ch.observe_algebra_element(claimed);
        let (proof, point) = prover.prove::<F, _>(&mut ch, n, 2, 0, claimed);

        let mut vch = challenger();
        let (vpoint, final_sum) = proof.verify(&mut vch, n, 2, 0).expect("verify");
        assert_eq!(vpoint.as_slice(), point.as_slice());

        let expected: F = originals
            .iter()
            .map(|(a, b)| {
                eval_mle(a.as_slice(), point.as_slice()) * eval_mle(b.as_slice(), point.as_slice())
            })
            .sum();
        assert_eq!(final_sum, expected);
    }

    /// Same for degree 3.
    #[test]
    fn deg3_roundtrip() {
        let mut rng = ChaCha8Rng::seed_from_u64(2);
        let n = 5;
        let triples = alloc::vec![(
            rand_poly(&mut rng, n),
            rand_poly(&mut rng, n),
            rand_poly(&mut rng, n),
        )];
        let originals = triples.clone();

        let mut prover = Deg3Prover::new(triples);
        let claimed = prover.claimed_sum();

        let mut ch = challenger();
        ch.observe_algebra_element(claimed);
        let (proof, point) = prover.prove::<F, _>(&mut ch, n, 3, 0, claimed);

        let mut vch = challenger();
        let (_, final_sum) = proof.verify(&mut vch, n, 3, 0).expect("verify");

        let (a, b, c) = &originals[0];
        let expected = eval_mle(a.as_slice(), point.as_slice())
            * eval_mle(b.as_slice(), point.as_slice())
            * eval_mle(c.as_slice(), point.as_slice());
        assert_eq!(final_sum, expected);
    }

    /// A tampered claimed sum changes the derived challenges, so the terminal
    /// sum can no longer be the honest product — which is the check the caller
    /// is responsible for.
    #[test]
    fn deg2_tampered_claim_breaks_the_terminal_check() {
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        let n = 4;
        let pairs = alloc::vec![(rand_poly(&mut rng, n), rand_poly(&mut rng, n))];
        let (a, b) = pairs[0].clone();

        let mut prover = Deg2Prover::new(pairs);
        let claimed = prover.claimed_sum();
        let mut ch = challenger();
        ch.observe_algebra_element(claimed);
        let (mut proof, _) = prover.prove::<F, _>(&mut ch, n, 2, 0, claimed);

        proof.claimed_sum += F::ONE;
        let mut vch = challenger();
        let (point, final_sum) = proof.verify(&mut vch, n, 2, 0).expect("shape is valid");
        let truth =
            eval_mle(a.as_slice(), point.as_slice()) * eval_mle(b.as_slice(), point.as_slice());
        assert_ne!(final_sum, truth);
    }

    /// Powers helper sanity, since the batching challenges lean on it.
    #[test]
    fn powers_are_powers() {
        let x = F::from_u8(5);
        assert_eq!(powers(x, 4), alloc::vec![F::ONE, x, x * x, x * x * x]);
    }
}
