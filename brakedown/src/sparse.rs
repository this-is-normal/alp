//! A minimal row-list sparse matrix, enough for the Brakedown expander blocks.
//!
//! Replaces the deleted `p3_matrix::sparse::CsrMatrix` + `mul_csr_dense`. Stored
//! as one `(column, value)` list per row — simple and clear rather than packed
//! CSR; the expander blocks are small and this is not a hot path in M2.

use alloc::vec;
use alloc::vec::Vec;

use p3_field::{Algebra, Field};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;
use rand::distr::{Distribution, StandardUniform};
use rand::{Rng, RngExt};

/// A sparse `height × width` matrix over `F`, stored as one entry-list per row.
#[derive(Clone, Debug)]
pub struct SparseMatrix<F> {
    width: usize,
    /// `rows[i]` = the nonzero `(column, value)` pairs of row `i`.
    rows: Vec<Vec<(usize, F)>>,
}

impl<F: Field> SparseMatrix<F> {
    /// Build from explicit rows. Panics if any column index is out of range.
    pub fn new(width: usize, rows: Vec<Vec<(usize, F)>>) -> Self {
        assert!(
            rows.iter()
                .all(|row| row.iter().all(|&(col, _)| col < width)),
            "sparse matrix column index out of range"
        );
        Self { width, rows }
    }

    /// Number of rows.
    pub fn height(&self) -> usize {
        self.rows.len()
    }

    /// Number of columns.
    pub const fn width(&self) -> usize {
        self.width
    }

    /// A random `height × width` matrix where every row has exactly `row_weight`
    /// nonzero entries, at distinct uniformly-random columns with uniform nonzero
    /// values. (The Brakedown expander blocks are fixed-row-weight sparse.)
    pub fn rand_fixed_row_weight<R: Rng>(
        rng: &mut R,
        height: usize,
        width: usize,
        row_weight: usize,
    ) -> Self
    where
        StandardUniform: Distribution<F>,
    {
        assert!(row_weight <= width, "row weight exceeds width");
        let rows = (0..height)
            .map(|_| {
                // Pick `row_weight` distinct columns.
                let mut cols: Vec<usize> = Vec::with_capacity(row_weight);
                while cols.len() < row_weight {
                    let c = (rng.random::<u64>() as usize) % width;
                    if !cols.contains(&c) {
                        cols.push(c);
                    }
                }
                // Assign each a uniform nonzero value.
                cols.into_iter()
                    .map(|c| {
                        let mut v: F = rng.random();
                        while v.is_zero() {
                            v = rng.random();
                        }
                        (c, v)
                    })
                    .collect()
            })
            .collect();
        Self { width, rows }
    }

    /// Compute `self · dense`, where `dense` has `self.width()` rows. The result
    /// has `self.height()` rows and `dense.width()` columns.
    ///
    /// Generic over the dense value type `A` (an `F`-algebra), so the same
    /// `F`-valued sparse matrix multiplies both base-field dense matrices
    /// (`A = F`, the committed rows) and extension-field ones (`A = EF`, the
    /// combined rows produced during opening). Each product `s * coeff` pairs a
    /// dense value `s: A` with a sparse coefficient `coeff: F` via `A`'s
    /// `Mul<F>`.
    pub fn mul_dense<A: Algebra<F> + Copy + Send + Sync>(
        &self,
        dense: &RowMajorMatrix<A>,
    ) -> RowMajorMatrix<A> {
        assert_eq!(
            dense.height(),
            self.width,
            "sparse·dense dimension mismatch"
        );
        let w = dense.width();
        let mut out = vec![A::ZERO; self.rows.len() * w];
        for (i, row) in self.rows.iter().enumerate() {
            let dst = &mut out[i * w..(i + 1) * w];
            for &(j, coeff) in row {
                let src = &dense.values[j * w..(j + 1) * w];
                for (d, &s) in dst.iter_mut().zip(src) {
                    // s : A (dense value), coeff : F — A is an F-algebra so `A * F` is defined.
                    *d += s * coeff;
                }
            }
        }
        RowMajorMatrix::new(out, w)
    }
}
