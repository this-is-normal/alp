//! **Quasar** — a field-agnostic multilinear polynomial commitment scheme with a
//! polylogarithmic verifier, from quasi-abelian codes
//! ([2026/1839](https://eprint.iacr.org/2026/1839)).
//!
//! Quasar keeps the fast commitment of a Brakedown-style tensor commitment over
//! a [QA code](p3_qa_code) — reshape, encode every row, Merkle-commit the
//! columns — and replaces Brakedown's `√N` proof and verifier with **code
//! switching** into [BaseFold](p3_basefold).
//!
//! # The opening protocol
//!
//! Two row combinations are formed: a proximity row `v = rᵗM` for random `r`,
//! and the evaluation row `p = z_Lᵗ M`. Both are QA-encoded and their `c` blocks
//! committed with BaseFold. Then:
//!
//! 1. **Selector sumcheck** — `t` sampled columns of the committed matrix are
//!    opened, and one degree-2 sumcheck over `log(c·N₂)` variables ties them to
//!    the committed blocks (both branches batched by one challenge).
//! 2. **Outer encoding sumcheck** — a degree-3 sumcheck over `n₂` variables
//!    certifies that every parity block really is `WHT(e_i ⋆ WHT(input))`. Its
//!    three factors are `H̃(γ,·)`, the batched *public* coefficients, and the
//!    encoding witness `WHT(input)`, which is **never committed**.
//! 3. **Inner encoding sumcheck** — a degree-2 sumcheck reduces the witness
//!    claim to an evaluation of the committed input block, using
//!    `WHT(v)~(r) = Σ_y H̃(r,y)·ṽ(y)`.
//! 4. **Claim reduction + two batched openings** — every claim so far, plus the
//!    evaluation claim `q₁(z_R) = y` (free: the code is systematic), is
//!    collapsed onto one point `ζ` by a degree-2 sumcheck, then discharged by
//!    one batched BaseFold opening of the online blocks and one of the
//!    preprocessed coefficients.
//!
//! The verifier evaluates `H̃` and the sparse selector itself in `O(log N)`
//! time; everything else it checks is either a sumcheck transcript or a
//! BaseFold opening, hence `O(λ log² N)` overall.
//!
//! # Instantiation
//!
//! Challenges live in the committed field itself (no extension), matching the
//! paper's Mersenne127 setting and this fork's ed25519 target. QA codes need a
//! large field for their distance bound — see [`p3_security::QaRegime`] — so
//! that is not a practical restriction: 31-bit fields are unusable for QA
//! regardless.

#![no_std]

extern crate alloc;

mod params;
mod pcs;
mod prescribed_pcs;
mod sumchecks;
mod util;

pub use params::QuasarParams;
pub use pcs::{QuasarError, QuasarPcs, QuasarProof, QuasarProverData, QuasarProverKey};
pub use prescribed_pcs::{QuasarMlError, QuasarMlPcs, QuasarMlProof, QuasarMlProverData};
pub use util::{eq_eval, eq_table, eval_mle, wht_mle_eval, wht_mle_table};
