//! One module per commitment-scheme backend, selected by cargo feature.
//!
//! Each backend exports the same two things — `CodePcs` (the multilinear PCS
//! `multi-stark` will use) and `build_pcs(log_height, width)` (one instance
//! sized for a trace of that shape). `config.rs` therefore needs no
//! per-backend code at all.
//!
//! | feature | scheme | proof / verifier |
//! |---|---|---|
//! | *(none)* | Brakedown tensor PCS | `O(√N)` |
//! | `ligero` | Reed–Solomon tensor PCS | `O(√N)`, high distance ⇒ few queries |
//! | `basefold` | BaseFold over a random foldable code | `O(λ log²N)` |
//! | `quasar` | Quasar: QA code switched into BaseFold | `O(λ log²N)`, much faster prover |

#[cfg(feature = "basefold")]
mod basefold;
#[cfg(all(not(feature = "basefold"), not(feature = "quasar"), feature = "ligero"))]
mod ligero;
#[cfg(feature = "quasar")]
mod quasar;
#[cfg(all(
    not(feature = "basefold"),
    not(feature = "quasar"),
    not(feature = "ligero")
))]
mod brakedown;

#[cfg(feature = "basefold")]
pub use basefold::{CodePcs, NAME, build_pcs};
#[cfg(all(not(feature = "basefold"), not(feature = "quasar"), feature = "ligero"))]
pub use ligero::{CodePcs, NAME, build_pcs};
#[cfg(feature = "quasar")]
pub use quasar::{CodePcs, NAME, build_pcs, resolved_params};
#[cfg(all(
    not(feature = "basefold"),
    not(feature = "quasar"),
    not(feature = "ligero")
))]
pub use brakedown::{CodePcs, NAME, build_pcs};

// `MultiStarkConfig::committed_table` needs no per-backend code: the
// `ProverData<Self>` projection normalises to the selected backend's concrete
// type, whose inherent `table()` then resolves.

/// Padding floor. Every backend here commits at the trace's natural arity, so
/// any positive arity is fine.
pub const FOLDING: usize = 1;
