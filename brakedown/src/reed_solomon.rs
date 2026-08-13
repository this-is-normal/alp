//! A naive, field-agnostic Reed–Solomon code — the base case of the Brakedown
//! recursion, and the "simple code" the tensor PCS (M3a) commits through
//! (an RS-tensor commitment is exactly Ligero).
//!
//! Encoding evaluates each column's message polynomial (coefficients = the
//! `k` message rows) at `n` distinct field points via a Vandermonde matrix.
//! `O(n·k)` per column, fine for the small instances, and needs only `n`
//! distinct field elements — no FFT, no two-adicity.

use alloc::vec;
use alloc::vec::Vec;

use p3_field::{Algebra, ExtensionField, Field};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;
use p3_maybe_rayon::prelude::*;

use crate::{LinearCode, LinearCodeExt};

/// A `[n, k]` Reed–Solomon code evaluated at the points `0, 1, …, n-1`.
#[derive(Clone, Debug)]
pub struct ReedSolomonCode<F> {
    /// Message dimension (polynomial degree bound).
    k: usize,
    /// The `n` distinct evaluation points.
    points: Vec<F>,
}

impl<F: Field> ReedSolomonCode<F> {
    /// A `[n, k]` RS code over the points `0..n`.
    ///
    /// Requires `n >= k` and that the field has at least `n` elements (so the
    /// points are distinct).
    pub fn new(k: usize, n: usize) -> Self {
        assert!(
            n >= k,
            "RS codeword length must be at least the message length"
        );
        let points = (0..n).map(F::from_usize).collect();
        Self { k, points }
    }
}

/// Vandermonde encoding generic over the message value type `A` (an `F`-algebra),
/// so the same `F`-linear code encodes both base-field messages (`A = F`, the
/// committed rows) and extension-field messages (`A = EF`, the combined rows).
///
/// Treats `msg` as `k` message rows × `w` independent columns and returns the
/// `n` × `w` codeword: `out[t][c] = Σ_i msg[i][c] · point_t^i`.
fn vandermonde_encode<F: Field, A: Algebra<F> + Copy + Send + Sync>(
    k: usize,
    points: &[F],
    msg: &RowMajorMatrix<A>,
) -> RowMajorMatrix<A> {
    assert_eq!(msg.height(), k, "message must have `k` rows");
    let w = msg.width();
    let n = points.len();
    let mut out = vec![A::ZERO; n * w];
    // Each evaluation point produces one independent output row (`out[t] = P(p_t)`),
    // so the outer loop parallelizes cleanly (rayon under the `parallel` feature,
    // serial otherwise). This is the dominant cost when RS is used as the
    // top-level Ligero code (`O(n·k)` per column).
    out.par_chunks_exact_mut(w)
        .zip(points.par_iter())
        .for_each(|(dst, &p)| {
            let mut pow = F::ONE;
            for i in 0..k {
                let src = &msg.values[i * w..(i + 1) * w];
                for (d, &s) in dst.iter_mut().zip(src) {
                    // s : A (message value), pow : F — A is an F-algebra so `A * F` is defined.
                    *d += s * pow;
                }
                pow *= p;
            }
        });
    RowMajorMatrix::new(out, w)
}

impl<F: Field> LinearCode<F> for ReedSolomonCode<F> {
    fn message_len(&self) -> usize {
        self.k
    }

    fn codeword_len(&self) -> usize {
        self.points.len()
    }

    fn encode(&self, message: &RowMajorMatrix<F>) -> RowMajorMatrix<F> {
        vandermonde_encode(self.k, &self.points, message)
    }
}

impl<F: Field, EF: ExtensionField<F>> LinearCodeExt<F, EF> for ReedSolomonCode<F> {
    fn encode_ext(&self, message: &RowMajorMatrix<EF>) -> RowMajorMatrix<EF> {
        vandermonde_encode(self.k, &self.points, message)
    }
}
