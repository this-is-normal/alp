//! A trivial, field-agnostic multilinear polynomial commitment scheme.
//!
//! This exists as an **honest-prover reference oracle** for the `multi-stark`
//! (SuperSpartan-style) proving path. Its entire reason to exist:
//!
//! 1. It requires only `Val: Field` — **no `TwoAdicField`**, no FFT — so it
//!    proves the multilinear proving path is genuinely field-agnostic.
//! 2. It is obviously correct, so it can be used to differential-test a real
//!    field-agnostic PCS (Brakedown) by comparing opened evaluations.
//!
//! ## What it is (and is NOT)
//!
//! `commit` produces a real, binding Merkle commitment over the committed
//! tables (via any [`Mmcs`]) and absorbs it into the transcript, exactly like a
//! real PCS. But `open_at` simply *computes* the requested multilinear
//! evaluations directly from the retained tables, and `verify_at` **trusts the
//! revealed evaluations** (it re-absorbs them for Fiat-Shamir consistency and
//! shape-checks them, but does not re-open the commitment to bind them to the
//! committed data).
//!
//! Therefore this PCS is sound **only against an honest prover**. It must not be
//! used to test opening-soundness (a malicious prover could open arbitrary
//! values). It is a correctness oracle, not a secure commitment. `multi-stark`'s
//! zerocheck still catches constraint/public-value violations independently of
//! the PCS; only *opening tampering* would go undetected here.

#![no_std]

extern crate alloc;

use alloc::vec::Vec;
use core::marker::PhantomData;

use p3_challenger::{CanObserve, CanSampleUniformBits, FieldChallenger, GrindingChallenger};
use p3_commit::{Mmcs, MultilinearPcs};
use p3_field::{ExtensionField, Field};
use p3_matrix::dense::RowMajorMatrix;
use p3_multilinear_util::point::Point;
use p3_sumcheck::layout::Table;
use p3_sumcheck::{OpeningBatch, OpeningEvals, OpeningProtocol, PrescribedPointPcs};
use serde::{Deserialize, Serialize};

/// A trivial reveal-everything multilinear PCS, generic over the base field
/// `Val`, challenge field `Challenge`, and a Merkle commitment scheme `M`.
#[derive(Clone)]
pub struct TrivialMultilinearPcs<Val, Challenge, M> {
    /// Merkle scheme used only to produce a real, observable commitment.
    mmcs: M,
    /// Stacked variable count reported by [`MultilinearPcs::num_vars`].
    /// Informational only; `multi-stark` does not consult it.
    num_vars: usize,
    _phantom: PhantomData<(Val, Challenge)>,
}

impl<Val, Challenge, M> TrivialMultilinearPcs<Val, Challenge, M> {
    /// Build a trivial PCS from a Merkle scheme and a (nominal) variable count.
    pub const fn new(mmcs: M, num_vars: usize) -> Self {
        Self {
            mmcs,
            num_vars,
            _phantom: PhantomData,
        }
    }
}

/// Opening-verification failures for [`TrivialMultilinearPcs`].
#[derive(Debug)]
pub enum TrivialPcsError {
    /// The proof carried a different number of opening batches than the protocol requested.
    OpeningCountMismatch { expected: usize, got: usize },
    /// An opening batch's per-side shape disagreed with the protocol request.
    OpeningShapeMismatch { batch: usize },
}

/// Opening proof: just the claimed per-batch evaluations (current then next).
#[derive(Clone, Serialize, Deserialize)]
pub struct TrivialProof<Challenge> {
    evals: Vec<OpeningBatch<Challenge>>,
}

impl<Val, Challenge, M, Challenger> MultilinearPcs<Challenge, Challenger>
    for TrivialMultilinearPcs<Val, Challenge, M>
where
    Val: Field,
    Challenge: ExtensionField<Val>,
    M: Mmcs<Val>,
    Challenger: FieldChallenger<Val> + CanObserve<M::Commitment>,
{
    type Val = Val;
    type Commitment = M::Commitment;
    type ProverData = Vec<Table<Val>>;
    type Proof = TrivialProof<Challenge>;
    type Error = TrivialPcsError;
    type Witness = Vec<Table<Val>>;
    type OpeningProtocol = OpeningProtocol;

    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn commit(
        &self,
        witness: Self::Witness,
        challenger: &mut Challenger,
    ) -> (Self::Commitment, Self::ProverData) {
        // Commit each table as one matrix (one row per column/poly, 2^arity wide).
        // The MMCS handles tables of differing heights in one commitment.
        let matrices: Vec<RowMajorMatrix<Val>> = witness
            .iter()
            .map(|t| {
                let width = 1usize << t.num_variables();
                let values: Vec<Val> = t.iter_polys().flatten().copied().collect();
                RowMajorMatrix::new(values, width)
            })
            .collect();
        let (commitment, _mmcs_data) = self.mmcs.commit(matrices);
        // Bind the commitment into the transcript; the multi-stark verifier
        // absorbs the same commitment before calling `verify_at`.
        challenger.observe(commitment.clone());
        (commitment, witness)
    }

    fn open(
        &self,
        _prover_data: Self::ProverData,
        _protocol: Self::OpeningProtocol,
        _challenger: &mut Challenger,
    ) -> Self::Proof {
        unimplemented!("TrivialMultilinearPcs supports only prescribed-point opening (open_at)")
    }

    fn verify(
        &self,
        _commitment: &Self::Commitment,
        _proof: &Self::Proof,
        _challenger: &mut Challenger,
        _protocol: Self::OpeningProtocol,
    ) -> Result<(), Self::Error> {
        unimplemented!("TrivialMultilinearPcs supports only prescribed-point verification (verify_at)")
    }
}

impl<Val, Challenge, M, Challenger> PrescribedPointPcs<Challenge, Challenger>
    for TrivialMultilinearPcs<Val, Challenge, M>
where
    Val: Field,
    Challenge: ExtensionField<Val>,
    M: Mmcs<Val>,
    Challenger: FieldChallenger<Val>
        + GrindingChallenger<Witness = Val>
        + CanSampleUniformBits<Val>
        + CanObserve<M::Commitment>,
{
    fn open_at(
        &self,
        prover_data: Self::ProverData,
        protocol: &OpeningProtocol,
        points: &[Point<Challenge>],
        challenger: &mut Challenger,
    ) -> Self::Proof {
        assert_eq!(
            protocol.num_openings(),
            points.len(),
            "one prescribed point per opening batch"
        );
        let mut evals = Vec::with_capacity(points.len());
        for ((table_idx, batch), point) in protocol.iter_openings().zip(points) {
            let table = &prover_data[table_idx];
            // Direct evaluations at the prescribed (table-local) point.
            let current: Vec<Challenge> = batch
                .current()
                .iter()
                .map(|&col| table.poly(col).eval_base(point))
                .collect();
            // Repeat-last successor ("next row") evaluations at the same point.
            let next: Vec<Challenge> = batch
                .next()
                .iter()
                .map(|&col| table.poly(col).eval_next_base(point))
                .collect();
            // Absorb current then next — the exact order `verify_at` replays.
            challenger.observe_algebra_slice(&current);
            challenger.observe_algebra_slice(&next);
            evals.push(OpeningBatch::new(current, next));
        }
        TrivialProof { evals }
    }

    fn verify_at(
        &self,
        _commitment: &Self::Commitment,
        proof: &Self::Proof,
        protocol: &OpeningProtocol,
        _points: &[Point<Challenge>],
        challenger: &mut Challenger,
    ) -> Result<Vec<OpeningEvals<Challenge>>, Self::Error> {
        assert_eq!(
            protocol.num_openings(),
            _points.len(),
            "one prescribed point per opening batch"
        );
        if proof.evals.len() != protocol.num_openings() {
            return Err(TrivialPcsError::OpeningCountMismatch {
                expected: protocol.num_openings(),
                got: proof.evals.len(),
            });
        }
        // Re-absorb the claimed evals in the same order the prover did, so the
        // transcript stays in lockstep. (Honest-oracle: the evals are trusted.)
        for (i, ((_table_idx, batch), claimed)) in
            protocol.iter_openings().zip(&proof.evals).enumerate()
        {
            if !batch.has_same_shape(claimed) {
                return Err(TrivialPcsError::OpeningShapeMismatch { batch: i });
            }
            challenger.observe_algebra_slice(claimed.current());
            challenger.observe_algebra_slice(claimed.next());
        }
        Ok(proof.evals.clone())
    }
}
