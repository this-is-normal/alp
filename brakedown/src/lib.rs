//! Brakedown linear code — the Spielman-based recursive linear-time code from
//! the Brakedown paper, rebuilt against current Plonky3 traits.
//!
//! This is the **encoder** (M2). It is field-agnostic (`Val: Field`, no FFT /
//! two-adicity) and is the linear code a Brakedown-style multilinear PCS
//! (M3/M4) will commit through.
//!
//! ## Encoding (batched, along the row axis)
//!
//! A message is a [`RowMajorMatrix`] with `message_len` rows and any number of
//! columns; each column is an independent length-`message_len` message, all
//! encoded together. The output has `codeword_len` rows, same width.
//!
//! The recursive Spielman construction ([`BrakedownCode`]):
//!
//! ```text
//!   encode(x):
//!     y = A · x            # A: sparse, fixed-row-weight random
//!     z = inner.encode(y)  # recurse; base case = Reed–Solomon
//!     v = B · z            # B: sparse
//!     return [ x ; z ; v ] # systematic: message x preserved on top
//! ```
//!
//! ## Scope (M2)
//!
//! Mechanics-first: correct, small-instance, explicit per-level parameters
//! (like the historical `brakedown!` macro). Production parameter generation
//! (distances/rates per the paper) and soundness accounting are deferred.

#![no_std]

extern crate alloc;

mod brakedown_code;
mod pcs;
mod prescribed_pcs;
mod reed_solomon;
mod sparse;

pub use brakedown_code::{BrakedownCode, BrakedownParams};
pub use pcs::{TensorPcs, TensorProof, TensorProverData, TensorPcsError};
pub use prescribed_pcs::{
    BrakedownMlPcs, BrakedownPcsError, BrakedownProof, BrakedownProverData,
};
pub use reed_solomon::ReedSolomonCode;
pub use sparse::SparseMatrix;

use p3_field::{ExtensionField, Field};
use p3_matrix::dense::RowMajorMatrix;

/// A linear code that encodes a batch of messages along the row axis.
///
/// Object-safe so a [`BrakedownCode`]'s inner code can be a `Box<dyn LinearCode>`.
pub trait LinearCode<F: Field> {
    /// Number of rows in an input message (the code's dimension `k`).
    fn message_len(&self) -> usize;

    /// Number of rows in an output codeword (the code's length `n`).
    fn codeword_len(&self) -> usize;

    /// Encode `message` (which must have [`message_len`](Self::message_len) rows)
    /// into a codeword with [`codeword_len`](Self::codeword_len) rows and the same width.
    fn encode(&self, message: &RowMajorMatrix<F>) -> RowMajorMatrix<F>;
}

/// A linear code that can also encode messages whose entries live in an
/// extension field `EF`. The code is `F`-linear, so the same generator applies
/// componentwise to `EF` messages. Parameterized by `EF` (rather than a generic
/// method) so it stays object-safe.
///
/// The tensor PCS needs this: committed rows are over `F`, but the combined
/// rows produced during opening live in the challenge field `EF`.
pub trait LinearCodeExt<F: Field, EF: ExtensionField<F>>: LinearCode<F> {
    /// Encode an `EF`-valued message with the same `F`-linear code as [`LinearCode::encode`].
    fn encode_ext(&self, message: &RowMajorMatrix<EF>) -> RowMajorMatrix<EF>;
}
