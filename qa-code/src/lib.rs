//! Quasi-abelian (QA) codes — the commit-side code of the Quasar PCS
//! ([2026/1839](https://eprint.iacr.org/2026/1839)), and on its own the code of
//! a Brakedown-style tensor PCS ("QAPCS").
//!
//! A rank-1 systematic QA code of index `c` encodes a length-`N` message
//! (`N = 2ⁿ`) as `c` blocks of length `N`, the first being the message itself.
//! Encoding is two Walsh–Hadamard transforms and a pointwise product per parity
//! block: **`2(c−1)N log N` additions and `(c−1)N` multiplications**, no FFT and
//! no two-adicity, so it runs over any odd-characteristic field — including the
//! ed25519 base field this fork targets.
//!
//! ```text
//!   Enc(v) = ( v ‖ WHT(e_2 ⋆ WHT(v)) ‖ … ‖ WHT(e_c ⋆ WHT(v)) ),   e_i = 2^{-n}·WHT(a_i)
//! ```
//!
//! # Layering
//!
//! - [`QaCode`] implements [`p3_brakedown::LinearCode`] and
//!   [`p3_brakedown::LinearCodeExt`], so it drops into
//!   [`p3_brakedown::TensorPcs`] / [`p3_brakedown::BrakedownMlPcs`] unchanged —
//!   that composition *is* QAPCS.
//! - The proven distance and the query count it implies live in
//!   [`p3_security::QaRegime`], not here: they are soundness parameters, not
//!   code mechanics.
//! - Quasar (`p3-quasar`) adds the code-switching layer on top, and reuses
//!   [`QaCode::encode_vec`] / [`QaCode::middle_vec`] / [`QaCode::scaled_coeffs`]
//!   for its encoding sumcheck.
//!
//! # Distance needs a big field
//!
//! The bound of [LLX⁺26b] is vacuous for 31-bit fields at any useful length
//! (≈0 relative distance), marginal at 64 bits, ≈0.43 at 127 bits and ≈0.46 at
//! 255 bits (rate ½). QA codes are a large-field construction.

#![no_std]

extern crate alloc;

mod code;
mod wht;

pub use code::{QaCode, generator_matrix};
pub use wht::{wht, wht_batched, wht_naive};
