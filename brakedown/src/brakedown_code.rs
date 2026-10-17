//! The recursive Spielman-based Brakedown code.

use alloc::boxed::Box;
use alloc::vec::Vec;

use p3_field::{Algebra, ExtensionField, Field};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;
use rand::Rng;
use rand::distr::{Distribution, StandardUniform};

use crate::{LinearCode, LinearCodeExt, ReedSolomonCode, SparseMatrix};

/// The Brakedown code: `encode(x) = [ x ; inner.encode(A·x) ; B·inner.encode(A·x) ]`.
///
/// Systematic — the message `x` is the top `message_len` rows of the codeword.
///
/// Generic over both the base field `F` (the sparse blocks `A`, `B` live over
/// `F`) and an extension field `EF`: the inner code is a
/// [`LinearCodeExt<F, EF>`], so the whole recursive code encodes both
/// base-field messages ([`LinearCode::encode`]) and extension-field messages
/// ([`LinearCodeExt::encode_ext`]) with the same `F`-linear generator.
pub struct BrakedownCode<F, EF> {
    a: SparseMatrix<F>,
    b: SparseMatrix<F>,
    inner: Box<dyn LinearCodeExt<F, EF>>,
}

impl<F: Field, EF: ExtensionField<F>> BrakedownCode<F, EF> {
    /// Assemble one Brakedown level from its two sparse blocks and an inner code.
    ///
    /// Dimensions must chain: `A` maps the message (`A.width()` rows) to the
    /// inner code's input (`A.height() == inner.message_len()`), and `B` maps
    /// the inner codeword (`B.width() == inner.codeword_len()`) to the tail block.
    pub fn new(
        a: SparseMatrix<F>,
        b: SparseMatrix<F>,
        inner: Box<dyn LinearCodeExt<F, EF>>,
    ) -> Self {
        assert_eq!(
            a.height(),
            inner.message_len(),
            "A.height() must equal inner.message_len() (A·x feeds the inner code)"
        );
        assert_eq!(
            b.width(),
            inner.codeword_len(),
            "B.width() must equal inner.codeword_len() (B·z consumes the inner codeword)"
        );
        Self { a, b, inner }
    }

    /// Generic encoder shared by [`encode`](LinearCode::encode) and
    /// [`encode_ext`](LinearCodeExt::encode_ext).
    ///
    /// The recursion `[ x ; z ; v ]` with `y = A·x`, `z = inner_encode(y)`,
    /// `v = B·z` is identical for both value types `A` (`= F` or `= EF`); only
    /// the inner encode differs (base-field vs. extension-field), passed in as a
    /// closure so the sparse-multiply / stack logic is written once.
    fn encode_with<T, InnerEnc>(&self, x: &RowMajorMatrix<T>, inner_encode: InnerEnc) -> RowMajorMatrix<T>
    where
        T: Algebra<F> + Copy + Send + Sync,
        InnerEnc: FnOnce(&RowMajorMatrix<T>) -> RowMajorMatrix<T>,
    {
        assert_eq!(
            x.height(),
            self.message_len(),
            "message must have `message_len` rows"
        );
        let w = x.width();

        let y = self.a.mul_dense(x); // A·x
        let z = inner_encode(&y); // inner codeword of y
        let v = self.b.mul_dense(&z); // B·z

        // Stack [ x ; z ; v ] vertically (all share width `w`, so row blocks concatenate).
        let mut values: Vec<T> = Vec::with_capacity((x.height() + z.height() + v.height()) * w);
        values.extend_from_slice(&x.values);
        values.extend_from_slice(&z.values);
        values.extend_from_slice(&v.values);
        RowMajorMatrix::new(values, w)
    }
}

impl<F: Field, EF: ExtensionField<F>> LinearCode<F> for BrakedownCode<F, EF> {
    fn message_len(&self) -> usize {
        self.a.width()
    }

    fn codeword_len(&self) -> usize {
        // [ x (message) ; z (inner codeword) ; v (B·z) ]
        self.message_len() + self.inner.codeword_len() + self.b.height()
    }

    fn encode(&self, x: &RowMajorMatrix<F>) -> RowMajorMatrix<F> {
        self.encode_with(x, |y| self.inner.encode(y))
    }
}

impl<F: Field, EF: ExtensionField<F>> LinearCodeExt<F, EF> for BrakedownCode<F, EF> {
    fn encode_ext(&self, x: &RowMajorMatrix<EF>) -> RowMajorMatrix<EF> {
        self.encode_with(x, |y| self.inner.encode_ext(y))
    }
}

/// Structural parameters for the recursive Brakedown construction.
///
/// # Caveat on the constants (M4 scope boundary)
///
/// These are **structurally-correct** Brakedown parameters — the recursion
/// shape, the fixed-row-weight sparse expander blocks, and the rate/dimension
/// relations all follow the Brakedown paper (Golovnev, Lee, Setty, Thaler,
/// Wahby, *"Brakedown: Linear-time and field-agnostic SNARKs for R1CS"*). The
/// degree (row-weight) and rate constants here are **reasonable defaults, not
/// formally-derived** distance parameters: a provable minimum relative distance
/// requires the paper's exact degree formulas (its `c_n`, `d_n` derived from the
/// target distance/false-positive rates and the expander mixing lemma). Using
/// them for a concrete security level is deferred — the code is *linear-time and
/// correct*, but its distance is assumed, not proven. See also the soundness note
/// in [`crate::TensorPcs`].
#[derive(Clone, Debug)]
pub struct BrakedownParams {
    /// Recursion base: message lengths `<= base_threshold` bottom out in a
    /// Reed–Solomon code instead of recursing.
    pub base_threshold: usize,
    /// Inner message length `m = round(alpha * message_len)`, as `alpha_num / alpha_den`.
    pub alpha_num: usize,
    pub alpha_den: usize,
    /// Codeword length `n = round(rate * message_len)`, as `rate_num / rate_den`.
    pub rate_num: usize,
    pub rate_den: usize,
    /// Fixed row weight of the first expander block `A`.
    pub c_a: usize,
    /// Fixed row weight of the second expander block `B`.
    pub c_b: usize,
}

impl Default for BrakedownParams {
    /// `alpha ≈ 0.238`, `rate ≈ 1.72`, base threshold `16`, row weights `7`/`10`.
    fn default() -> Self {
        Self {
            base_threshold: 16,
            alpha_num: 238,
            alpha_den: 1000,
            rate_num: 172,
            rate_den: 100,
            c_a: 7,
            c_b: 10,
        }
    }
}

/// `round(x * num / den)` in integer arithmetic (no `std` float rounding).
fn round_ratio(x: usize, num: usize, den: usize) -> usize {
    (x * num + den / 2) / den
}

impl<F, EF> BrakedownCode<F, EF>
where
    F: Field,
    EF: ExtensionField<F>,
    StandardUniform: Distribution<F>,
{
    /// Build the recursive Brakedown code for any `message_len`, using the
    /// [default parameters](BrakedownParams::default).
    ///
    /// The top-level `message_len` will typically be a power of two (the tensor
    /// PCS needs `cols = message_len = 2^j`); intermediate inner lengths need
    /// not be. See [`BrakedownParams`] for the important caveat that the
    /// distance constants are structural defaults, not formally derived.
    pub fn from_message_len(message_len: usize, rng: &mut impl Rng) -> Self {
        Self::from_message_len_with_params(message_len, &BrakedownParams::default(), rng)
    }

    /// Build the recursive Brakedown code for any `message_len` with explicit
    /// [`BrakedownParams`].
    ///
    /// Structure (per the Brakedown paper):
    /// * **base** (`message_len <= params.base_threshold`): a Reed–Solomon code
    ///   `RS(message_len, round(rate * message_len))`, which is a genuine
    ///   `LinearCodeExt`.
    /// * **recurse**: inner message length `m = round(alpha * message_len)`;
    ///   `A` is a random `m × message_len` fixed-row-weight (`c_a`) sparse block;
    ///   the inner code is built recursively on `m`; `B` is a random
    ///   `b_rows × inner_n` fixed-row-weight (`c_b`) block with `b_rows` chosen so
    ///   the total `codeword_len = message_len + inner_n + b_rows ≈ round(rate * message_len)`.
    pub fn from_message_len_with_params(
        message_len: usize,
        params: &BrakedownParams,
        rng: &mut impl Rng,
    ) -> Self {
        assert!(message_len > 0, "message_len must be positive");

        // Inner message length m = round(alpha * message_len); clamp to a valid,
        // strictly-smaller length so the recursion always makes progress.
        let m = round_ratio(message_len, params.alpha_num, params.alpha_den)
            .max(1)
            .min(message_len - 1);

        // A: m × message_len, fixed row weight c_a (clamped to the width).
        let c_a = params.c_a.min(message_len).max(1);
        let a = SparseMatrix::<F>::rand_fixed_row_weight(rng, m, message_len, c_a);

        // Inner code on the length-m message.
        let inner: Box<dyn LinearCodeExt<F, EF>> = if m <= params.base_threshold {
            let n_base = round_ratio(m, params.rate_num, params.rate_den).max(m);
            Box::new(ReedSolomonCode::<F>::new(m, n_base))
        } else {
            Box::new(Self::from_message_len_with_params(m, params, rng))
        };
        let inner_n = inner.codeword_len();

        // B: b_rows × inner_n, with b_rows set so the total codeword length hits
        // the target rate. Guaranteed positive for the default rate (rate·ml >
        // ml + inner_n), but saturate for safety at tiny sizes.
        let target = round_ratio(message_len, params.rate_num, params.rate_den);
        let b_rows = target.saturating_sub(message_len + inner_n).max(1);
        let c_b = params.c_b.min(inner_n).max(1);
        let b = SparseMatrix::<F>::rand_fixed_row_weight(rng, b_rows, inner_n, c_b);

        Self::new(a, b, inner)
    }
}
