//! The rank-1 systematic quasi-abelian code.

use alloc::vec;
use alloc::vec::Vec;

use p3_brakedown::{LinearCode, LinearCodeExt};
use p3_field::{Algebra, ExtensionField, Field};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;
use p3_maybe_rayon::prelude::*;
use rand::{Rng, RngExt};
use rand::distr::{Distribution, StandardUniform};
use serde::{Deserialize, Serialize};

use crate::wht::wht_batched;

/// A rank-1 systematic QA code of index `c` (= inverse rate) and message
/// length `2^log_message_len`.
///
/// # The code
///
/// Identify a message `v ∈ F^N` (`N = 2ⁿ`) with an element of the group algebra
/// `R = F[X₁,…,X_n]/(X₁²−1, …, X_n²−1)`. For public `a₁,…,a_c ∈ R` the codeword
/// is `(v·a₁, …, v·a_c)`; systematic means `a₁ = 1`, so the first block is the
/// message. Multiplication in `R` diagonalizes under the Walsh–Hadamard
/// transform, giving the encoder
///
/// ```text
///   v'    = WHT(v)
///   u_i   = WHT( (2^{-n}·a_i') ⋆ v' )      for i = 2..c,   a_i' = WHT(a_i)
///   Enc(v) = ( v ‖ u_2 ‖ … ‖ u_c )
/// ```
///
/// so encoding costs `2(c−1)·N log N` additions and `(c−1)·N` multiplications —
/// no FFT, no two-adicity, and no field inverse (the `2^{-n}` is folded into the
/// stored coefficients).
///
/// The stored `coeffs[i]` *is* `2^{-n}·a_{i+2}'`; sampling those uniformly is
/// equivalent to sampling `a_i` uniformly, since scaling by a fixed nonzero
/// constant is a bijection.
///
/// # Distance
///
/// The relative distance of a *random* QA code is a random variable; use
/// [`p3_security::QaRegime`] for the proven bound and the query count it
/// implies. It is only useful over large fields: vacuous for 31-bit fields,
/// ≈0.43 at 127 bits, ≈0.46 at 255 bits (rate ½).
///
/// # Field requirement
///
/// Odd characteristic (`2` must be invertible). Characteristic 2 collapses the
/// transform, and [`QaCode::new_random`] rejects it.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QaCode<F> {
    log_message_len: usize,
    /// `c − 1` coefficient vectors of length `2^log_message_len`, one per parity
    /// block; already scaled by `2^{-log_message_len}`.
    coeffs: Vec<Vec<F>>,
}

impl<F: Field> QaCode<F> {
    /// A random QA code of index `index` (≥ 2) and message length
    /// `2^log_message_len`.
    pub fn new_random<R: Rng>(log_message_len: usize, index: usize, rng: &mut R) -> Self
    where
        StandardUniform: Distribution<F>,
    {
        assert!(index >= 2, "QA index (inverse rate) must be at least 2");
        assert!(
            (F::ONE + F::ONE) != F::ZERO,
            "QA codes need odd characteristic (2 must be invertible)"
        );
        let n = 1usize << log_message_len;
        let coeffs = (0..index - 1)
            .map(|_| (0..n).map(|_| rng.random::<F>()).collect())
            .collect();
        Self {
            log_message_len,
            coeffs,
        }
    }

    /// Build from explicit *scaled* coefficient vectors `2^{-n}·a_i'`, one per
    /// parity block.
    pub fn from_scaled_coeffs(log_message_len: usize, coeffs: Vec<Vec<F>>) -> Self {
        let n = 1usize << log_message_len;
        assert!(!coeffs.is_empty(), "QA index must be at least 2");
        assert!(
            coeffs.iter().all(|c| c.len() == n),
            "each coefficient vector must have the message length"
        );
        Self {
            log_message_len,
            coeffs,
        }
    }

    /// `log₂` of the message length.
    pub const fn log_message_len(&self) -> usize {
        self.log_message_len
    }

    /// The QA index `c` = inverse rate.
    pub const fn index(&self) -> usize {
        self.coeffs.len() + 1
    }

    /// The scaled coefficient vectors `2^{-n}·a_i'`, `i = 2..c`.
    ///
    /// Quasar commits to these during setup and opens them at the challenge
    /// points of the encoding sumcheck, so they are public parameters.
    pub fn scaled_coeffs(&self) -> &[Vec<F>] {
        &self.coeffs
    }

    /// Encode a single message vector, returning the `c·N` codeword.
    ///
    /// The one-row case of [`LinearCode::encode`]; the sumcheck layers of
    /// Quasar work per row, so they use this directly.
    pub fn encode_vec<A: Algebra<F> + Copy + Send + Sync>(&self, message: &[A]) -> Vec<A> {
        assert_eq!(message.len(), self.message_len(), "wrong message length");
        let mut out = Vec::with_capacity(self.codeword_len());
        out.extend_from_slice(message);

        let mut middle = message.to_vec();
        wht_batched(&mut middle, 1);
        self.append_parity_blocks(&middle, 1, &mut out);
        out
    }

    /// `WHT(v)` for a message vector — the intermediate the encoding sumcheck
    /// treats as a witness (it is never committed).
    pub fn middle_vec<A: Algebra<F> + Copy + Send + Sync>(&self, message: &[A]) -> Vec<A> {
        assert_eq!(message.len(), self.message_len(), "wrong message length");
        let mut middle = message.to_vec();
        wht_batched(&mut middle, 1);
        middle
    }

    /// Shared encoder body: `[ message ; parity blocks ]`, generic over the
    /// value type so one implementation serves base-field and extension-field
    /// messages.
    fn encode_matrix<A: Algebra<F> + Copy + Send + Sync>(
        &self,
        message: &RowMajorMatrix<A>,
    ) -> RowMajorMatrix<A> {
        assert_eq!(
            message.height(),
            self.message_len(),
            "message must have `message_len` rows"
        );
        let width = message.width();

        let mut middle = message.values.clone();
        wht_batched(&mut middle, width);

        let mut values = Vec::with_capacity(self.codeword_len() * width);
        values.extend_from_slice(&message.values);
        self.append_parity_blocks(&middle, width, &mut values);

        debug_assert_eq!(values.len(), self.codeword_len() * width);
        RowMajorMatrix::new(values, width)
    }

    /// `u_i = WHT(coeffs[i] ⋆ middle)` for each parity block, appended to `out`.
    ///
    /// `middle` is `WHT(message)` laid out as `message_len × width` row-major;
    /// the coefficient for row `r` multiplies that whole row.
    fn append_parity_blocks<A: Algebra<F> + Copy + Send + Sync>(
        &self,
        middle: &[A],
        width: usize,
        out: &mut Vec<A>,
    ) {
        let n = self.message_len();
        for coeff in &self.coeffs {
            let start = out.len();
            out.resize(start + n * width, A::ZERO);
            let block = &mut out[start..];
            block
                .par_chunks_exact_mut(width)
                .zip(middle.par_chunks_exact(width))
                .zip(coeff.par_iter())
                .for_each(|((dst, src), &e)| {
                    for (d, &s) in dst.iter_mut().zip(src) {
                        // s : A (message value), e : F — A is an F-algebra.
                        *d = s * e;
                    }
                });
            wht_batched(block, width);
        }
    }
}

impl<F: Field> LinearCode<F> for QaCode<F> {
    fn message_len(&self) -> usize {
        1 << self.log_message_len
    }

    fn codeword_len(&self) -> usize {
        self.index() << self.log_message_len
    }

    fn encode(&self, message: &RowMajorMatrix<F>) -> RowMajorMatrix<F> {
        self.encode_matrix(message)
    }
}

impl<F: Field, EF: ExtensionField<F>> LinearCodeExt<F, EF> for QaCode<F> {
    fn encode_ext(&self, message: &RowMajorMatrix<EF>) -> RowMajorMatrix<EF> {
        self.encode_matrix(message)
    }
}

/// The generator matrix `(I | A_2 | … | A_c)` of this code, for tests.
///
/// `A_i` is the dyadic matrix of `a_i`, i.e. `A_i[r][s] = a_i[r ⊕ s]`, and
/// `a_i = WHT(2^{-n}·a_i')` — recovering `a_i` from the stored scaled spectrum
/// is exactly one more WHT. `2^{codeword bits}` entries; small instances only.
pub fn generator_matrix<F: Field>(code: &QaCode<F>) -> RowMajorMatrix<F> {
    let n = code.message_len();
    let c = code.index();
    let mut values = vec![F::ZERO; n * c * n];

    // Systematic block.
    for r in 0..n {
        values[r * (c * n) + r] = F::ONE;
    }

    for (i, scaled) in code.scaled_coeffs().iter().enumerate() {
        let mut a = scaled.clone();
        wht_batched(&mut a, 1);
        for r in 0..n {
            for s in 0..n {
                values[r * (c * n) + (i + 1) * n + s] = a[r ^ s];
            }
        }
    }
    RowMajorMatrix::new(values, c * n)
}
