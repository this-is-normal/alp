//! A **foldable random linear code** (Basefold, ed25519-native — no FFT).
//!
//! # Construction (rate `1/c`)
//! Base code `C_0 = [c, 1]`: a message `[x]` encodes to `x · g0` for a fixed random `g0 ∈ (F*)^c`.
//! Level `i` from level `i-1`, with a random "twist" `t_i ∈ (F*)^{n_{i-1}}` (a diagonal):
//! ```text
//!   Enc_i(m_L ‖ m_R) = ( Enc_{i-1}(m_L) + t_i · Enc_{i-1}(m_R) ,
//!                        Enc_{i-1}(m_L) − t_i · Enc_{i-1}(m_R) )
//! ```
//! so a length-`2^i` message maps to a length-`c·2^i` codeword (rate `1/c` at every level).
//!
//! # Folding (MLE convention)
//! One butterfly layer inverts with the round challenge `α`. Writing a codeword `w = (a ‖ b)`
//! (halves of length `n_{i-1}`), recover `l = Enc(m_L)`, `r = Enc(m_R)` then take the MLE
//! restriction of the folded variable to `α`:
//! ```text
//!   l = (a + b)/2,  r = (a − b)/(2·t_i),  fold(w, α) = l + α·(r − l)   [ = (1−α)·l + α·r ]
//! ```
//! The load-bearing identity (Basefold soundness rests on it) is **commutativity** with the
//! multilinear restriction of the folded variable:
//! ```text
//!   fold( Enc_i(m), α ) = Enc_{i-1}( (1−α)·m_L + α·m_R )
//! ```
//! so folding the codeword by `α` equals encoding the message with its top variable bound to `α`,
//! matching what a sumcheck round does. M1's test checks this. The verifier's per-query check
//! recomputes `fold` at a single index using `t_i[j]`.

use alloc::vec::Vec;

use p3_field::{Field, batch_multiplicative_inverse};
use p3_maybe_rayon::prelude::*;
use rand::distr::{Distribution, StandardUniform};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;

/// Recursion level below which `encode_level` stops forking: the subproblem is
/// then small enough that task overhead outweighs the parallelism.
const PAR_ENCODE_MIN_LEVEL: usize = 8;

/// Odd stride (golden-ratio constant) decorrelating the per-element seeds.
const SEED_STRIDE: u64 = 0x9E37_79B9_7F4A_7C15;

/// Sample the nonzero code element at global index `idx`, **addressably**.
///
/// The construction used to draw every element from one sequential ChaCha
/// stream, which meant a verifier could only reach `twist[level][j]` by
/// replaying every earlier draw — i.e. by materialising the whole `O(N)` code
/// to check `O(q·log N)` positions. Seeding a stream per element instead makes
/// each position reachable on its own, which is what
/// [`RandomFoldableCode::new_verifier`] needs. Rejection stays inside the
/// element's own stream, so it does not disturb addressability.
fn sample_nonzero_at<F: Field>(seed: u64, idx: usize) -> F
where
    StandardUniform: Distribution<F>,
{
    let mut rng =
        ChaCha8Rng::seed_from_u64(seed.wrapping_add((idx as u64).wrapping_mul(SEED_STRIDE)));
    loop {
        let x: F = rng.random();
        if !x.is_zero() {
            return x;
        }
    }
}

/// Global index of `twist[level][0]`. `g0` occupies `0..inv_rate`, and level
/// `i` the `inv_rate·2^{i−1}` slots starting at `inv_rate·2^{i−1}`.
const fn level_base(inv_rate: usize, level: usize) -> usize {
    inv_rate << (level - 1)
}

/// The interface a Basefold code must provide.
pub trait FoldableCode<F: Field> {
    /// Number of foldable rounds = log2 of the message length.
    fn log_len(&self) -> usize;
    /// Length of a full codeword (`c · 2^{log_len}`).
    fn codeword_len(&self) -> usize;
    /// Encode a full message of length `2^{log_len}` into a codeword of length `codeword_len`.
    fn encode(&self, message: &[F]) -> Vec<F>;
    /// Fold a codeword of message-log-size `level` by one round with challenge `alpha`,
    /// halving its length. `level` must be in `1..=log_len`.
    fn fold(&self, codeword: &[F], level: usize, alpha: F) -> Vec<F>;
    /// One entry of the twist vector `t_level`, the only twist access the
    /// verifier's fold-consistency check needs.
    fn diag_at(&self, level: usize, j: usize) -> F;
}

/// A rate-`1/inv_rate` foldable random linear code over `F`, built deterministically from a seed
/// (ChaCha, so the construction is byte-portable across host/zkVM).
#[derive(Clone, Debug)]
pub struct RandomFoldableCode<F> {
    log_len: usize,
    inv_rate: usize,
    /// The seed every element is derived from, kept so a lazy instance can
    /// reach any position without storing it.
    seed: u64,
    /// Base `[inv_rate, 1]` generator; all entries nonzero.
    g0: Vec<F>,
    /// `twist[i]` (for `i in 1..=log_len`) has length `inv_rate · 2^{i-1}`; all entries nonzero.
    /// `twist[0]` is unused (kept for 1-based indexing).
    /// `None` on a lazy verifier instance — see [`Self::new_verifier`].
    twist: Option<Vec<Vec<F>>>,
    /// Elementwise inverses of `twist`, likewise absent when lazy.
    tinv: Option<Vec<Vec<F>>>,
    /// Lazy sampler, present exactly when `twist` is absent.
    ///
    /// A plain function pointer rather than a `StandardUniform: Distribution<F>`
    /// bound on the struct: the bound would have to be repeated on every impl
    /// that touches a code, and from there on `BasefoldPcs`, `QuasarPcs` and
    /// each `MultilinearPcs` impl above them. The constructor already carries
    /// the bound, so it can hand the sampler over as data.
    sample: Option<fn(u64, usize) -> F>,
    /// `1/2` in `F`.
    inv2: F,
}

impl<F: Field> RandomFoldableCode<F>
where
    StandardUniform: Distribution<F>,
{
    /// Build the full code for `2^log_len`-element messages at rate `1/inv_rate`,
    /// seeded by `seed`. This is the **prover's** constructor: it materialises
    /// every twist. See [`Self::new_verifier`] for the other side.
    pub fn new(log_len: usize, inv_rate: usize, seed: u64) -> Self {
        assert!(inv_rate >= 2, "inverse rate must be >= 2");
        let g0: Vec<F> = (0..inv_rate).map(|i| sample_nonzero_at(seed, i)).collect();

        let mut twist = Vec::with_capacity(log_len + 1);
        let mut tinv = Vec::with_capacity(log_len + 1);
        twist.push(Vec::new()); // index 0 unused
        tinv.push(Vec::new());
        for i in 1..=log_len {
            let n_prev = inv_rate << (i - 1); // c · 2^{i-1}
            // Addressable sampling, so every draw is independent and the level
            // can be filled in parallel.
            let base = level_base(inv_rate, i);
            let t: Vec<F> = (0..n_prev)
                .into_par_iter()
                .map(|j| sample_nonzero_at(seed, base + j))
                .collect();
            // Montgomery's trick: ONE inversion plus ~3n multiplications for the
            // whole level, instead of n independent inversions. Over ed25519 an
            // inversion is a^(p-2) — 252 squarings and 77 multiplications — so
            // inverting `c·2^log_len` twists one at a time dominated setup
            // entirely (355 s at log_len = 19; this makes it ~1 s).
            tinv.push(batch_multiplicative_inverse(&t));
            twist.push(t);
        }

        let inv2 = F::TWO.try_inverse().expect("2 is invertible in an odd field");
        Self { log_len, inv_rate, seed, g0, twist: Some(twist), tinv: Some(tinv), sample: None, inv2 }
    }

    /// The **verifier's** constructor: the same code, holding only the seed and
    /// the base generator.
    ///
    /// A Basefold verifier reads the twists at `O(q·log N)` queried positions,
    /// but building them all is `O(N)` — sampling `inv_rate·(2^log_len − 1)`
    /// field elements. On a host that sets up once this is amortised away; in a
    /// recursion circuit, which sets up per proof, it *is* the cost, and it
    /// grows linearly while the verifier proper grows logarithmically.
    ///
    /// [`Self::diag_at`] then derives one twist on demand from
    /// [`sample_nonzero_at`], and `fold_pair` inverts that single element
    /// rather than a whole level.
    pub fn new_verifier(log_len: usize, inv_rate: usize, seed: u64) -> Self {
        assert!(inv_rate >= 2, "inverse rate must be >= 2");
        let g0: Vec<F> = (0..inv_rate).map(|i| sample_nonzero_at(seed, i)).collect();
        let inv2 = F::TWO.try_inverse().expect("2 is invertible in an odd field");
        Self {
            log_len,
            inv_rate,
            seed,
            g0,
            twist: None,
            tinv: None,
            sample: Some(sample_nonzero_at::<F>),
            inv2,
        }
    }

}

impl<F: Field> RandomFoldableCode<F> {
    /// `1/t_level[j]`: an array read when the code is materialised, one
    /// sample-and-invert when it is lazy.
    fn tinv_at(&self, level: usize, j: usize) -> F {
        match &self.tinv {
            Some(tinv) => tinv[level][j],
            None => self
                .diag_at(level, j)
                .try_inverse()
                .expect("twists are sampled nonzero"),
        }
    }

    /// The materialised twists. Encoding and whole-codeword folding are prover
    /// operations; a lazy verifier instance has neither.
    ///
    /// # Panics
    ///
    /// If this instance came from [`Self::new_verifier`].
    fn materialised_twist(&self) -> &[Vec<F>] {
        self.twist
            .as_deref()
            .expect("this is a lazy (verifier-only) code: it cannot encode")
    }

    /// As [`Self::materialised_twist`], for the inverses.
    fn materialised_tinv(&self) -> &[Vec<F>] {
        self.tinv
            .as_deref()
            .expect("this is a lazy (verifier-only) code: it cannot fold a whole codeword")
    }

    /// Encode a message of length `2^level` (`level in 0..=log_len`). Public for tests; the trait's
    /// `encode` calls this at the top level.
    pub fn encode_level(&self, level: usize, message: &[F]) -> Vec<F> {
        debug_assert_eq!(message.len(), 1usize << level, "message length must be 2^level");
        if level == 0 {
            return self.g0.iter().map(|&g| message[0] * g).collect();
        }
        let half = message.len() / 2;
        // The two halves are independent, so fork them until the subproblems
        // are too small to be worth a task. This recursion IS the encode, and
        // BaseFold commits a full-length codeword with it, so it dominates the
        // commitment for the standalone scheme.
        let (l, r) = if level >= PAR_ENCODE_MIN_LEVEL {
            p3_maybe_rayon::prelude::join(
                || self.encode_level(level - 1, &message[..half]),
                || self.encode_level(level - 1, &message[half..]),
            )
        } else {
            (
                self.encode_level(level - 1, &message[..half]),
                self.encode_level(level - 1, &message[half..]),
            )
        };
        let t = &self.materialised_twist()[level];
        debug_assert_eq!(l.len(), t.len());
        let n_prev = l.len();

        // out = [ l + t·r ; l − t·r ], both halves written in one parallel pass.
        let mut out = F::zero_vec(2 * n_prev);
        let (top, bottom) = out.split_at_mut(n_prev);
        top.par_iter_mut()
            .zip(bottom.par_iter_mut())
            .zip(l.par_iter().zip(r.par_iter()).zip(t.par_iter()))
            .for_each(|((hi, lo), ((&lj, &rj), &tj))| {
                let tr = tj * rj;
                *hi = lj + tr;
                *lo = lj - tr;
            });
        out
    }

    /// Fold a single butterfly pair `(a, b) = (O[j], O[j+n/2])` at `level` and index `j`,
    /// matching `fold(codeword, level, alpha)[j]`. Used by the verifier's per-query check.
    pub fn fold_pair(&self, a: F, b: F, level: usize, j: usize, alpha: F) -> F {
        let l = (a + b) * self.inv2;
        let r = (a - b) * self.inv2 * self.tinv_at(level, j);
        l + alpha * (r - l)
    }

    /// Base-code encoding of a single value `v`: `[v·g0[0], …, v·g0[c−1]]` (the final oracle `O_d`).
    pub fn base_encode(&self, v: F) -> Vec<F> {
        self.g0.iter().map(|&g| v * g).collect()
    }

    /// Fold a message (not a codeword) by one round (MLE restrict of the top variable to `alpha`):
    /// `[m_L ‖ m_R] -> (1−α)·m_L + α·m_R`. For tests.
    pub fn fold_message(message: &[F], alpha: F) -> Vec<F> {
        let half = message.len() / 2;
        (0..half).map(|j| message[j] + alpha * (message[half + j] - message[j])).collect()
    }
}

impl<F: Field> FoldableCode<F> for RandomFoldableCode<F> {
    fn log_len(&self) -> usize {
        self.log_len
    }

    fn codeword_len(&self) -> usize {
        self.inv_rate << self.log_len
    }

    fn encode(&self, message: &[F]) -> Vec<F> {
        assert_eq!(message.len(), 1usize << self.log_len, "message length must be 2^log_len");
        self.encode_level(self.log_len, message)
    }

    fn fold(&self, codeword: &[F], level: usize, alpha: F) -> Vec<F> {
        assert!((1..=self.log_len).contains(&level), "level out of range");
        let expected = self.inv_rate << level;
        assert_eq!(codeword.len(), expected, "codeword length must be inv_rate·2^level");
        let n_prev = codeword.len() / 2;
        let (a, b) = codeword.split_at(n_prev);
        let tinv = &self.materialised_tinv()[level];
        // Every output position is independent.
        a.par_iter()
            .zip(b.par_iter())
            .zip(tinv.par_iter())
            .map(|((&aj, &bj), &tij)| {
                let l = (aj + bj) * self.inv2;
                let r = (aj - bj) * self.inv2 * tij;
                l + alpha * (r - l) // MLE restrict: (1−α)·l + α·r
            })
            .collect()
    }

    fn diag_at(&self, level: usize, j: usize) -> F {
        match (&self.twist, self.sample) {
            (Some(twist), _) => twist[level][j],
            (None, Some(sample)) => sample(self.seed, level_base(self.inv_rate, level) + j),
            (None, None) => unreachable!("a code has either its twists or a sampler"),
        }
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use p3_ed25519::Ed25519;
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    use super::*;

    type F = Ed25519;

    fn rand_vec(rng: &mut ChaCha8Rng, n: usize) -> Vec<F> {
        (0..n).map(|_| rng.random::<F>()).collect()
    }

    /// The load-bearing identity: folding the CODEWORD equals encoding the folded MESSAGE.
    /// Checked at every level, for several random messages/challenges.
    #[test]
    fn fold_commutes_with_encode() {
        let d = 6;
        let code = RandomFoldableCode::<F>::new(d, 4, 0xBA5E_F01D);
        let mut rng = ChaCha8Rng::seed_from_u64(1);
        for _ in 0..8 {
            // Start from a full message and fold down level by level.
            let mut msg = rand_vec(&mut rng, 1 << d);
            let mut cw = code.encode_level(d, &msg);
            for level in (1..=d).rev() {
                let alpha: F = rng.random();
                let folded_cw = code.fold(&cw, level, alpha);
                let folded_msg = RandomFoldableCode::<F>::fold_message(&msg, alpha);
                let reencoded = code.encode_level(level - 1, &folded_msg);
                assert_eq!(folded_cw, reencoded, "fold∘encode ≠ encode∘fold at level {level}");
                msg = folded_msg;
                cw = folded_cw;
            }
            // After d folds we are at the base code C_0: a single message element -> c symbols.
            assert_eq!(cw.len(), code.inv_rate);
        }
    }

    /// `encode` is a deterministic linear map: same seed + message ⇒ identical codeword,
    /// and it is F-linear (Enc(a·x + b·y) = a·Enc(x) + b·Enc(y)).
    #[test]
    fn encode_is_deterministic_and_linear() {
        let d = 5;
        let code = RandomFoldableCode::<F>::new(d, 4, 42);
        let code2 = RandomFoldableCode::<F>::new(d, 4, 42);
        let mut rng = ChaCha8Rng::seed_from_u64(7);
        let x = rand_vec(&mut rng, 1 << d);
        let y = rand_vec(&mut rng, 1 << d);
        let a: F = rng.random();
        let b: F = rng.random();

        assert_eq!(code.encode(&x), code2.encode(&x), "same seed must give same code");

        let cx = code.encode(&x);
        let cy = code.encode(&y);
        let comb: Vec<F> = x.iter().zip(&y).map(|(&xi, &yi)| a * xi + b * yi).collect();
        let lhs = code.encode(&comb);
        let rhs: Vec<F> = cx.iter().zip(&cy).map(|(&u, &v)| a * u + b * v).collect();
        assert_eq!(lhs, rhs, "encode must be F-linear");
    }

    /// The lazy verifier code is the *same* code: every twist it derives on
    /// demand equals the materialised one, and so does every `fold_pair`.
    #[test]
    fn lazy_matches_full() {
        let (log_len, inv_rate, seed) = (6, 4, 0xBA5E_F01D);
        let full = RandomFoldableCode::<F>::new(log_len, inv_rate, seed);
        let lazy = RandomFoldableCode::<F>::new_verifier(log_len, inv_rate, seed);

        assert_eq!(full.g0, lazy.g0, "base generator");
        assert_eq!(full.codeword_len(), lazy.codeword_len());

        let mut rng = ChaCha8Rng::seed_from_u64(7);
        for level in 1..=log_len {
            for j in 0..(inv_rate << (level - 1)) {
                assert_eq!(
                    full.diag_at(level, j),
                    lazy.diag_at(level, j),
                    "twist mismatch at level {level} index {j}"
                );
            }
            // The per-query check the verifier actually runs.
            let (a, b, alpha) = (rng.random::<F>(), rng.random::<F>(), rng.random::<F>());
            assert_eq!(
                full.fold_pair(a, b, level, 0, alpha),
                lazy.fold_pair(a, b, level, 0, alpha),
                "fold_pair mismatch at level {level}"
            );
        }
    }
}