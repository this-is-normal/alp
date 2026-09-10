//! Basefold multilinear PCS over a **foldable random linear code** — an ed25519-native
//! (no-FFT) alternative to the tensor/Brakedown PCS with a polylog verifier.
//!
//! Build status (see `../BASEFOLD_PCS_PLAN.md`):
//! - **M1 (this crate, now):** the foldable code — [`RandomFoldableCode`] with `encode` and
//!   `fold`, satisfying the folding-commutativity relation that Basefold's soundness rests on.
//! - M2: the core commit/open/verify PCS.
//! - M3: `MultilinearPcs` + `PrescribedPointPcs` impls (drop-in for `BrakedownMlPcs`).

#![no_std]

extern crate alloc;

mod grinding;

mod foldable_code;
mod merkle_pcs;
mod pcs;
mod prescribed_pcs;

pub use foldable_code::{FoldableCode, RandomFoldableCode};
pub use grinding::grind_parallel;
pub use merkle_pcs::{BasefoldBatchProverData, BasefoldPcs, BasefoldProof, BasefoldProverData};
pub use pcs::{BasefoldCore, BasefoldCoreProof, eval_multilinear};
pub use prescribed_pcs::{BasefoldMlError, BasefoldMlPcs, BasefoldMlProof, BasefoldMlProverData};
