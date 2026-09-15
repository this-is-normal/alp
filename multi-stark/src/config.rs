//! Configuration tying together the commitment scheme, challenge field, and transcript.

use alloc::vec::Vec;

use p3_commit::MultilinearPcs;
use p3_field::{ExtensionField, Field};
pub use p3_multilinear_util::point::Point;
use p3_sumcheck::layout::Table;

/// The wiring a multilinear AIR proof depends on.
///
/// One implementation fixes a commitment scheme, a challenge field, and a transcript type.
///
/// The prover and verifier are generic over any implementation.
/// Current- and next-row evaluations of one table's preprocessed columns.
///
/// `next` carries **every** column; the verifier projects it onto the columns
/// the AIR actually reads a successor of.
#[derive(Clone, Debug)]
pub struct PreprocessedEvals<EF> {
    /// Current-row value of each column, in column order.
    pub current: Vec<EF>,
    /// Next-row value of each column, in column order.
    pub next: Vec<EF>,
}

pub trait MultiStarkConfig {
    /// Base field the trace and committed columns live in.
    type Val: Field;

    /// Extension field that challenges are drawn from.
    ///
    /// Large enough that a random challenge collides with negligible probability.
    type Challenge: ExtensionField<Self::Val>;

    /// Fiat-Shamir transcript used to derive challenges and absorb commitments.
    type Challenger;

    /// Multilinear commitment scheme over the base field.
    type Pcs: MultilinearPcs<Self::Challenge, Self::Challenger, Val = Self::Val>;

    /// Borrow the commitment scheme for the main trace.
    fn pcs(&self) -> &Self::Pcs;

    /// Borrow the commitment scheme for the preprocessed trace.
    ///
    /// The preprocessed trace stacks a different column count than the main trace.
    /// Its stacked polynomial therefore has a different arity.
    /// That arity needs its own scheme instance.
    ///
    /// Called only during setup, and only when the AIR declares preprocessed columns.
    /// Configurations for AIRs without a preprocessed trace need not override it.
    ///
    /// # Panics
    ///
    /// The default panics.
    /// A configuration proving an AIR with preprocessed columns must override it.
    fn preprocessed_pcs(&self) -> &Self::Pcs {
        unimplemented!("this configuration does not provide a preprocessed commitment scheme")
    }

    /// Smallest table arity the commitment scheme accepts without padding.
    ///
    /// A table below this floor is zero-padded before commitment, which breaks the successor view:
    /// - padding moves the repeated boundary row into the pad,
    /// - so the repeat-last successor view no longer reads itself,
    /// - which disagrees with the zerocheck successor convention.
    ///
    /// The committed prover therefore requires the trace arity to meet this floor.
    /// The opening point and the committed table then share one frame.
    ///
    /// For a WHIR-backed config this is the first-round folding factor.
    fn min_num_variables(&self) -> usize;

    /// Pack committed trace tables into the commitment scheme's witness form.
    ///
    /// The witness representation is private to the commitment scheme:
    /// - it folds in the slot layout and preprocessing depth that scheme expects,
    /// - the generic commitment interface exposes neither,
    /// - so the columns-to-witness bridge lives here, where those details are known.
    ///
    /// # Arguments
    ///
    /// - `tables`: source tables whose rows are trace-column polynomials.
    fn build_witness(
        &self,
        tables: Vec<Table<Self::Val>>,
    ) -> <Self::Pcs as MultilinearPcs<Self::Challenge, Self::Challenger>>::Witness;

    /// Borrows one committed trace table retained between commit and open.
    fn committed_table<'a>(
        &self,
        prover_data: &'a ProverData<Self>,
        table_index: usize,
    ) -> &'a Table<Self::Val>;
    /// Whether the verifier evaluates the preprocessed traces itself instead of
    /// opening a commitment to them.
    ///
    /// A preprocessed trace is public and fixed by the AIR, so a verifier that
    /// can evaluate it directly is in a *stronger* position than one trusting
    /// an opening — and it costs the protocol a commitment, an opening, and the
    /// proof bytes that carry it. Only worth it when the evaluation is cheap:
    /// `EdDsaAir`'s preprocessed trace is one 256-row pattern tiled per
    /// signature, so the work is constant in the batch size while the committed
    /// polynomial grows with it.
    ///
    /// Default `false` — commit and open as usual.
    fn evaluates_preprocessed_directly(&self) -> bool {
        false
    }

    /// Evaluate every preprocessed table at its point, in setup order over the
    /// AIRs that declare preprocessed columns.
    ///
    /// Called only when [`Self::evaluates_preprocessed_directly`] is true; the
    /// default panics so a config cannot enable the capability without
    /// supplying it.
    fn preprocessed_evals_at(
        &self,
        _points: &[Point<Self::Challenge>],
    ) -> Vec<PreprocessedEvals<Self::Challenge>> {
        unimplemented!(
            "config enables evaluates_preprocessed_directly but does not implement \
             preprocessed_evals_at"
        )
    }

}

/// Commitment scheme selected by a configuration.
pub type Pcs<C> = <C as MultiStarkConfig>::Pcs;

/// Succinct commitment produced by a configuration's commitment scheme.
pub type Commitment<C> = <Pcs<C> as MultilinearPcs<
    <C as MultiStarkConfig>::Challenge,
    <C as MultiStarkConfig>::Challenger,
>>::Commitment;

/// Prover-only data retained between commit and open.
pub type ProverData<C> = <Pcs<C> as MultilinearPcs<
    <C as MultiStarkConfig>::Challenge,
    <C as MultiStarkConfig>::Challenger,
>>::ProverData;

/// Opening proof produced by a configuration's commitment scheme.
pub type PcsProof<C> = <Pcs<C> as MultilinearPcs<
    <C as MultiStarkConfig>::Challenge,
    <C as MultiStarkConfig>::Challenger,
>>::Proof;

/// Opening-verification error produced by a configuration's commitment scheme.
pub type PcsError<C> = <Pcs<C> as MultilinearPcs<
    <C as MultiStarkConfig>::Challenge,
    <C as MultiStarkConfig>::Challenger,
>>::Error;
