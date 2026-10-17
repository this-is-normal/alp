//! Multi-stark backend for the EdDSA AIR, over the ed25519 base field.
//!
//! `EdDsaAir` declares a **preprocessed trace** (the fixed-base doubling ladder +
//! row selectors), so the config provides a *second* commitment scheme via
//! [`MultiStarkConfig::preprocessed_pcs`]; the preprocessed trace stacks a
//! different column count than the main trace, giving it a different stacked
//! arity, so it needs its own PCS instance sized for that arity.
//!
//! The field/hash/transcript stack lives in [`crate::hash`] and the commitment
//! scheme in [`crate::backend`] — one module per PCS, selected by cargo feature
//! (`ligero`, `basefold`, `quasar`, or none for Brakedown). Nothing in this file
//! is backend-specific.

use p3_multi_stark::config::{MultiStarkConfig, Point, PreprocessedEvals};
use p3_ed25519_air::eddsa::preprocessed_evals_at;
use p3_sumcheck::layout::{Layout, Table, Witness};

use crate::backend::{CodePcs, FOLDING, build_pcs};
use crate::hash::{EF, L, MyChallenger};

pub use crate::hash::{EF as ChallengeField, F, MyChallenger as Challenger, challenger};
// Only the ed25519-hash stack has hash muls to subtract; under `bb-hash` the
// hash is BabyBear, so there is nothing to calibrate.
#[cfg(all(feature = "count-ops", not(feature = "bb-hash")))]
pub use crate::hash::muls_per_perm;

/// A config for an AIR with a preprocessed trace: it carries two PCS instances,
/// one per stacked arity (main trace vs preprocessed trace).
pub struct EdDsaConfig {
    main_pcs: CodePcs,
    preprocessed_pcs: CodePcs,
    /// Scalar bit length, needed to rebuild the preprocessed pattern the
    /// verifier evaluates instead of opening.
    num_bits: usize,
}

impl MultiStarkConfig for EdDsaConfig {
    type Val = F;
    type Challenge = EF;
    type Challenger = MyChallenger;
    type Pcs = CodePcs;

    fn pcs(&self) -> &CodePcs {
        &self.main_pcs
    }

    fn preprocessed_pcs(&self) -> &CodePcs {
        &self.preprocessed_pcs
    }

    fn min_num_variables(&self) -> usize {
        FOLDING
    }

    fn build_witness(&self, tables: Vec<Table<F>>) -> Witness<F> {
        L::new_witness(tables, FOLDING)
    }

    /// The preprocessed trace is evaluated, not committed.
    ///
    /// `EdDsaAir`'s preprocessed columns are one 256-row pattern tiled once per
    /// signature, so their multilinear extension collapses to the pattern's and
    /// the verifier can evaluate it in ~5,600 multiplications — constant in the
    /// batch size, against a committed polynomial of 2^16 evaluations at 16
    /// signatures and 2^22 at 1024. It also removes ~42% of the proof.
    fn evaluates_preprocessed_directly(&self) -> bool {
        true
    }

    fn preprocessed_evals_at(
        &self,
        points: &[Point<EF>],
    ) -> Vec<PreprocessedEvals<EF>> {
        points
            .iter()
            .map(|point| {
                let (current, next) = preprocessed_evals_at(
                    crate::curve::d(),
                    crate::curve::base_b(),
                    self.num_bits,
                    point.as_slice(),
                );
                PreprocessedEvals { current, next }
            })
            .collect()
    }

    fn committed_table<'a>(
        &self,
        prover_data: &'a p3_multi_stark::config::ProverData<Self>,
        table_index: usize,
    ) -> &'a Table<F> {
        prover_data.table(table_index)
    }
}

impl EdDsaConfig {
    /// Both PCS instances, main first. A recursion host needs them to read out
    /// public setup material (Quasar's coefficient commitments, say) that the
    /// guest would otherwise have to recompute.
    pub const fn pcs_pair(&self) -> (&CodePcs, &CodePcs) {
        (&self.main_pcs, &self.preprocessed_pcs)
    }
}

/// Default target soundness in bits. Backends that take a caller-supplied
/// target (Brakedown, Ligero) derive their column-query count from it; the
/// others fix their own.
pub const TARGET_SECURITY_BITS: usize = 100;

/// Build the EdDSA config: both traces are `2^log_height` rows tall; the main
/// trace is `main_width` columns and the preprocessed trace `preprocessed_width`.
pub fn config_for_eddsa(
    num_bits: usize,
    log_height: usize,
    main_width: usize,
    preprocessed_width: usize,
) -> EdDsaConfig {
    config_for_eddsa_with_security(
        num_bits,
        log_height,
        main_width,
        preprocessed_width,
        TARGET_SECURITY_BITS,
    )
}

/// As [`config_for_eddsa`], but with an explicit soundness target.
///
/// Production uses [`TARGET_SECURITY_BITS`]. Small traces — tests, mostly —
/// may need a lower one: secure Brakedown at 100 bits needs enough columns to
/// open `t` queries, which a tiny `N` cannot supply, and `build_pcs` panics
/// rather than prove insecurely.
pub fn config_for_eddsa_with_security(
    num_bits: usize,
    log_height: usize,
    main_width: usize,
    preprocessed_width: usize,
    security_bits: usize,
) -> EdDsaConfig {
    EdDsaConfig {
        main_pcs: build_pcs(log_height, main_width, security_bits),
        preprocessed_pcs: build_pcs(log_height, preprocessed_width, security_bits),
        num_bits,
    }
}
