//! M2 (core): the Basefold commit / prove-evaluation / verify engine, validating the
//! **algebraic coupling** between the `eq(z,·)·f` sumcheck and the code fold chain.
//!
//! This core deliberately omits Merkle authentication and uses a minimal placeholder transcript:
//! it sends the intermediate oracles in full and checks the *entire* fold chain, so it proves the
//! protocol's algebra is right (honest opening verifies; a wrong value or a tampered oracle is
//! rejected). Merkle commitments + `q_bf`-query succinctness + the real Fiat–Shamir challenger are
//! layered on next (still M2), then the `MultilinearPcs`/`PrescribedPointPcs` impls (M3).

use alloc::vec;
use alloc::vec::Vec;

use p3_field::Field;
use p3_maybe_rayon::prelude::*;

use crate::foldable_code::{FoldableCode, RandomFoldableCode};

/// Minimal deterministic transcript (placeholder for a real challenger — M3). Both prover and
/// verifier drive it identically, which is all the M2-core correctness test needs.
struct Transcript<F> {
    state: F,
}
impl<F: Field> Transcript<F> {
    fn new() -> Self {
        Self { state: F::from_u64(0x9E37_79B9_7F4A_7C15) }
    }
    fn absorb(&mut self, xs: &[F]) {
        let c = F::from_u64(1_000_003);
        for &x in xs {
            self.state = self.state * c + x;
        }
    }
    fn challenge(&mut self) -> F {
        let c = F::from_u64(1_000_003);
        self.state = self.state * c + F::from_u64(0xC0FF_EE);
        self.state
    }
}

/// Full multilinear evaluation `f(z)` over the Boolean hypercube, MSB-first: index `x`'s bit at
/// position `d-1-k` is variable `k` (paired with `z[k]`).
pub fn eval_multilinear<F: Field>(evals: &[F], z: &[F]) -> F {
    let d = z.len();
    assert_eq!(evals.len(), 1 << d);
    eq_table(z).iter().zip(evals).map(|(&e, &f)| e * f).sum()
}

/// `eq_table(z)[x] = Π_k (x_bit_k ? z[k] : 1−z[k])`, MSB-first indexing (variable `k` at bit `d-1-k`).
pub(crate) fn eq_table<F: Field>(z: &[F]) -> Vec<F> {
    // Build by doubling, MSB-first: each variable in turn becomes the new
    // *least* significant bit, so entries for a variable land adjacent
    // (`next[2i]`, `next[2i+1]`) rather than half a table apart. Getting this
    // backwards silently reverses the variable order.
    let mut table = vec![F::ONE];
    for &zk in z {
        let mut next = F::zero_vec(2 * table.len());
        next.par_chunks_exact_mut(2)
            .zip(table.par_iter())
            .for_each(|(pair, &t)| {
                pair[0] = t * (F::ONE - zk);
                pair[1] = t * zk;
            });
        table = next;
    }
    table
}

/// MLE-restrict the top (MSB) variable to `r`: `[f_L ‖ f_R] -> (1−r)·f_L + r·f_R`.
pub(crate) fn bind_msb<F: Field>(f: &[F], r: F) -> Vec<F> {
    let half = f.len() / 2;
    let (lo, hi) = f.split_at(half);
    lo.par_iter()
        .zip(hi.par_iter())
        .map(|(&a, &b)| a + r * (b - a))
        .collect()
}

/// Evaluate a degree-2 univariate given by its values at 0,1,2, at point `r` (Lagrange).
pub(crate) fn interp_deg2<F: Field>(h: [F; 3], r: F) -> F {
    let inv2 = F::TWO.try_inverse().expect("2 invertible");
    let l0 = (r - F::ONE) * (r - F::TWO) * inv2;
    let l1 = r * (r - F::TWO);
    let l2 = r * (r - F::ONE) * inv2;
    h[0] * l0 - h[1] * l1 + h[2] * l2
}

/// A Basefold evaluation proof (M2-core: carries full intermediate oracles, no Merkle yet).
#[derive(Clone, Debug)]
pub struct BasefoldCoreProof<F> {
    /// Per round `i`: the sumcheck univariate `h_i` at points 0,1,2.
    pub sumcheck: Vec<[F; 3]>,
    /// Intermediate oracles `O_1 .. O_{d-1}` (folded codewords). `O_0` is the commitment,
    /// `O_d = v·g0` is derived from `final_value`.
    pub oracles: Vec<Vec<F>>,
    /// The fully folded message value `v = f(r)`.
    pub final_value: F,
}

/// The Basefold PCS core over one field (challenge field = base field; alp's eddsa config uses
/// `EF = Ed25519 = F`, so this covers the real use case).
pub struct BasefoldCore<F> {
    pub code: RandomFoldableCode<F>,
}

impl<F: Field> BasefoldCore<F> {
    pub fn new(code: RandomFoldableCode<F>) -> Self {
        Self { code }
    }

    fn num_vars(&self) -> usize {
        self.code.log_len()
    }

    /// Commit: encode the evaluations into the base oracle `O_0`.
    pub fn commit(&self, evals: &[F]) -> Vec<F> {
        assert_eq!(evals.len(), 1 << self.num_vars());
        self.code.encode(evals)
    }

    /// Prove `f(z) = y`. `o0` is the committed base oracle from [`commit`].
    pub fn prove(&self, evals: &[F], o0: &[F], z: &[F]) -> BasefoldCoreProof<F> {
        let d = self.num_vars();
        assert_eq!(z.len(), d);
        let mut tr = Transcript::<F>::new();
        tr.absorb(&[F::from_u64(o0.iter().fold(0u64, |a, _| a).wrapping_add(o0.len() as u64))]);

        // Product-sumcheck of g(x) = eq(z,x)·f(x): maintain BOTH multilinears and bind both by r_i.
        let mut f = evals.to_vec();
        let mut e = eq_table(z); // eq(z,·) over the hypercube
        let mut cw = o0.to_vec();
        let mut oracles: Vec<Vec<F>> = Vec::new(); // O_1..O_{d-1}
        let mut sumcheck: Vec<[F; 3]> = Vec::with_capacity(d);

        for i in 0..d {
            let half = f.len() / 2;
            let (fl, fr) = f.split_at(half);
            let (el, er) = e.split_at(half);
            // h_i(X) = Σ_{x'} E_X(x')·f_X(x'), degree 2; values at X = 0,1,2.
            let h0: F = el.iter().zip(fl).map(|(&a, &b)| a * b).sum();
            let h1: F = er.iter().zip(fr).map(|(&a, &b)| a * b).sum();
            let h2: F = (0..half)
                .map(|x| (F::TWO * er[x] - el[x]) * (F::TWO * fr[x] - fl[x]))
                .sum();
            let h = [h0, h1, h2];
            tr.absorb(&h);
            let r = tr.challenge();
            sumcheck.push(h);

            f = bind_msb(&f, r);
            e = bind_msb(&e, r);
            cw = self.code.fold(&cw, d - i, r);
            if i + 1 < d {
                oracles.push(cw.clone());
            }
        }

        let v = f[0]; // f(r)
        BasefoldCoreProof { sumcheck, oracles, final_value: v }
    }

    /// Verify `f(z) = y` against the base oracle `o0`.
    pub fn verify(&self, o0: &[F], z: &[F], y: F, proof: &BasefoldCoreProof<F>) -> Result<(), &'static str> {
        let d = self.num_vars();
        if proof.sumcheck.len() != d {
            return Err("wrong sumcheck length");
        }
        if proof.oracles.len() != d.saturating_sub(1) {
            return Err("wrong oracle count");
        }
        let mut tr = Transcript::<F>::new();
        tr.absorb(&[F::from_u64(o0.iter().fold(0u64, |a, _| a).wrapping_add(o0.len() as u64))]);

        // Re-derive challenges and check the sumcheck chain.
        let mut claim = y;
        let mut rs = Vec::with_capacity(d);
        for h in &proof.sumcheck {
            if h[0] + h[1] != claim {
                return Err("sumcheck round check failed");
            }
            tr.absorb(h);
            let r = tr.challenge();
            claim = interp_deg2(*h, r);
            rs.push(r);
        }
        // Final sumcheck claim must equal eq(z,r)·v.
        let eqzr = eval_eq(z, &rs);
        if claim != eqzr * proof.final_value {
            return Err("final sumcheck ≠ eq(z,r)·v");
        }

        // Fold chain: O_0 →(r_0) O_1 →…→ O_{d-1} →(r_{d-1}) O_d = v·g0.
        let mut cur = o0.to_vec();
        for i in 0..d {
            let folded = self.code.fold(&cur, d - i, rs[i]);
            let next = if i + 1 < d {
                proof.oracles[i].clone()
            } else {
                // O_d = base-code encoding of the single value v = v·g0.
                self.code.encode_level(0, &[proof.final_value])
            };
            if folded != next {
                return Err("fold chain mismatch");
            }
            cur = next;
        }
        Ok(())
    }
}

/// `eq(z, r) = Π_k eq(z[k], r[k])`.
pub(crate) fn eval_eq<F: Field>(z: &[F], r: &[F]) -> F {
    z.iter().zip(r).map(|(&zk, &rk)| zk * rk + (F::ONE - zk) * (F::ONE - rk)).product()
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use p3_ed25519::Ed25519;
    use p3_field::PrimeCharacteristicRing;
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    use super::*;

    type F = Ed25519;

    fn setup(d: usize) -> (BasefoldCore<F>, Vec<F>, Vec<F>, F, ChaCha8Rng) {
        let code = RandomFoldableCode::<F>::new(d, 4, 0xBA5E_F01D);
        let core = BasefoldCore::new(code);
        let mut rng = ChaCha8Rng::seed_from_u64(9);
        let evals: Vec<F> = (0..(1 << d)).map(|_| rng.random::<F>()).collect();
        let z: Vec<F> = (0..d).map(|_| rng.random::<F>()).collect();
        let y = eval_multilinear(&evals, &z);
        (core, evals, z, y, rng)
    }

    #[test]
    fn honest_open_verifies() {
        for d in 1..=8 {
            let (core, evals, z, y, _) = setup(d);
            let o0 = core.commit(&evals);
            let proof = core.prove(&evals, &o0, &z);
            // final_value is f(r) at the sumcheck challenges, NOT f(z); verify checks the
            // reduction f(z) = eq(z,r)·f(r) internally.
            core.verify(&o0, &z, y, &proof).unwrap_or_else(|e| panic!("honest reject at d={d}: {e}"));
        }
    }

    #[test]
    fn wrong_value_rejected() {
        let (core, evals, z, y, _) = setup(6);
        let o0 = core.commit(&evals);
        let proof = core.prove(&evals, &o0, &z);
        assert!(core.verify(&o0, &z, y + F::ONE, &proof).is_err(), "wrong y must be rejected");
    }

    #[test]
    fn tampered_oracle_rejected() {
        let (core, evals, z, y, _) = setup(6);
        let o0 = core.commit(&evals);
        let mut proof = core.prove(&evals, &o0, &z);
        // Corrupt one entry of an intermediate oracle -> fold chain must break.
        proof.oracles[0][0] += F::ONE;
        assert!(core.verify(&o0, &z, y, &proof).is_err(), "tampered oracle must be rejected");
    }
}
