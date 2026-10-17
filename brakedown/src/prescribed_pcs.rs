//! A multilinear PCS adapter (M3b) that makes the standalone tensor (Ligero)
//! commitment usable by the multivariate STARK (`p3-multi-stark`).
//!
//! # What this bridges
//!
//! [`TensorPcs`](crate::TensorPcs) proves a *single* point evaluation
//! `P(z) = v` of one committed multilinear polynomial. `multi-stark`'s opening
//! phase instead needs to prove *many* per-column claims (current and
//! repeat-last successor views, across several tables) on the *stacked*
//! committed polynomial. This adapter reduces the latter to the former with the
//! reusable stacked-layout sumcheck from `p3-sumcheck`, exactly the way the WHIR
//! adapter does — only the final `P(r) = v` step is discharged by the tensor
//! opening instead of WHIR's FRI.
//!
//! # The reduction (mirrors `whir/src/pcs/adapter.rs`)
//!
//! * **commit**: stack the tables into one polynomial `P` (via the layout
//!   [`Witness`]), tensor-commit `P`'s evaluations, and retain the layout for
//!   table access during opening.
//! * **open_at**: record every prescribed-point claim against a
//!   [`layout::Verifier`], batch them under a sampled `alpha` into one
//!   [`Constraint`] with weight polynomial `W`, run the reusable sumcheck prover
//!   to fold **all** stacked variables down to a point `r`, then tensor-open `P`
//!   at `r`. The proof bundles the recorded evals, the sumcheck transcript data,
//!   and the tensor opening proof.
//! * **verify_at**: rebuild the identical `Constraint` from the proof's evals,
//!   replay the sumcheck to recover `r` and the reduced claim, evaluate `W(r)`,
//!   tensor-verify the opening to recover `v = P(r)`, and accept iff
//!   `claimed_eval == W(r) · v`.
//!
//! Unlike WHIR there are **no out-of-domain (OOD) virtual evaluations**: the
//! tensor commitment carries its own proximity argument, so the reduction needs
//! only the honest opening claims.
//!
//! # Fiat-Shamir alignment
//!
//! Prover and verifier drive the *same* transcript operations in the same
//! order: claim recording absorbs the evals (both sides run the shared
//! [`layout::Verifier::add_claim_at`]), `alpha` is sampled once, the sumcheck
//! round loop observes/samples symmetrically, and the tensor opening's internal
//! `eta` / column-query challenges happen after the sumcheck on both sides. The
//! commitment is absorbed once — by the tensor `commit` on the prover side and
//! by the outer `multi-stark` verifier before `verify_at` — so `verify_at` (and
//! the tensor [`verify_opening`](TensorPcs::verify_opening) it calls) must not
//! re-absorb it.

use alloc::vec::Vec;
use core::marker::PhantomData;

use p3_challenger::{
    CanObserve, CanSampleBits, CanSampleUniformBits, FieldChallenger, GrindingChallenger,
};
use p3_commit::{Mmcs, MultilinearPcs};
use p3_field::{ExtensionField, Field};
use p3_multilinear_util::point::Point;
use p3_multilinear_util::poly::Poly;
use p3_sumcheck::layout::{Layout, Table, Verifier, Witness};
use p3_sumcheck::product_polynomial::ProductPolynomial;
use p3_sumcheck::strategy::{SumcheckProver, VariableOrder};
use p3_sumcheck::{
    OpeningBatch, OpeningEvals, OpeningProtocol, PrescribedPointPcs, SumcheckData, SumcheckError,
};
use serde::{Deserialize, Serialize};

use crate::pcs::{TensorProof, TensorProverData};
use crate::{LinearCodeExt, TensorPcs, TensorPcsError};

/// A Brakedown/Ligero-backed multilinear PCS for the multi-stark opening path.
///
/// Wraps a standalone [`TensorPcs`] and pairs it with a stacked-layout binding
/// mode `L` (e.g. `PrefixProver`) that supplies the reusable sumcheck reduction.
#[derive(Clone)]
pub struct BrakedownMlPcs<Val, EF, C, M, L> {
    /// The underlying single-point tensor commitment scheme.
    tensor: TensorPcs<Val, EF, C, M>,
    /// Marker tying this PCS to a stacked-layout binding mode; carries no state.
    _marker: PhantomData<L>,
}

impl<Val, EF, C, M, L> BrakedownMlPcs<Val, EF, C, M, L> {
    /// Wrap a tensor commitment scheme into the multilinear-PCS adapter.
    ///
    /// The tensor scheme's row/column split fixes the stacked variable count it
    /// can commit: `2^log_rows · code.message_len() == 2^num_variables`.
    pub const fn new(tensor: TensorPcs<Val, EF, C, M>) -> Self {
        Self {
            tensor,
            _marker: PhantomData,
        }
    }
}

/// Prover-side handoff between commit and open.
///
/// Holds the committed layout (for table access) and the tensor prover data
/// (for the final `P(r) = v` opening). Both describe the *same* stacked
/// polynomial.
pub struct BrakedownProverData<Val, EF, M, L>
where
    Val: Send + Sync + Clone,
    M: Mmcs<Val>,
{
    /// Stacked layout retaining the source tables and per-table placements.
    layout: L,
    /// Stacked committed polynomial evaluations, reused to build the sumcheck.
    stacked: Vec<Val>,
    /// Tensor prover data behind the commitment; consumed by the tensor opening.
    tensor: TensorProverData<Val, M>,
    /// Marker tying the data to its challenge field; carries no runtime state.
    _marker: PhantomData<EF>,
}

// Hand-written so the `Clone` bound lands on the tensor prover data (whose
// Merkle prover data associated type must be `Clone`) rather than on `M`.
impl<Val, EF, M, L> Clone for BrakedownProverData<Val, EF, M, L>
where
    Val: Send + Sync + Clone,
    M: Mmcs<Val>,
    L: Clone,
    TensorProverData<Val, M>: Clone,
{
    fn clone(&self) -> Self {
        Self {
            layout: self.layout.clone(),
            stacked: self.stacked.clone(),
            tensor: self.tensor.clone(),
            _marker: PhantomData,
        }
    }
}

impl<Val, EF, M, L> BrakedownProverData<Val, EF, M, L>
where
    Val: Field,
    EF: ExtensionField<Val>,
    M: Mmcs<Val>,
    L: Layout<Val, EF>,
{
    /// Returns source table `id` retained by the committed layout.
    pub fn table(&self, id: usize) -> &Table<Val> {
        self.layout.table(id)
    }
}

/// A prescribed-point opening proof.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(
    serialize = "Val: Serialize, EF: Serialize, M::Proof: Serialize",
    deserialize = "Val: serde::de::DeserializeOwned, EF: serde::de::DeserializeOwned, M::Proof: serde::de::DeserializeOwned"
))]
pub struct BrakedownProof<Val, EF, M>
where
    Val: Send + Sync + Clone,
    M: Mmcs<Val>,
{
    /// Per-batch claimed evaluations (current then next), in opening order.
    evals: Vec<OpeningEvals<EF>>,
    /// Sumcheck round transcript folding the batched claim to a point `r`.
    sumcheck: SumcheckData<Val, EF>,
    /// Tensor opening proof for `P(r) = v` at the folded point.
    tensor: TensorProof<Val, EF, M>,
}

/// Opening-verification failures for [`BrakedownMlPcs`].
#[derive(Debug)]
pub enum BrakedownPcsError<E> {
    /// The tensor `P(r) = v` opening failed to verify.
    Tensor(TensorPcsError<E>),
    /// The sumcheck round transcript was malformed (e.g. wrong round count).
    Sumcheck(SumcheckError),
    /// A claim batch's shape disagreed with the protocol request.
    OpeningShape(SumcheckError),
    /// The proof carried a different number of opening batches than requested.
    OpeningCountMismatch { expected: usize, got: usize },
    /// The reduced claim did not match `W(r) · v`.
    FinalCheckMismatch,
}

impl<EF, Val, C, M, L, Challenger> MultilinearPcs<EF, Challenger>
    for BrakedownMlPcs<Val, EF, C, M, L>
where
    Val: Field,
    EF: ExtensionField<Val>,
    C: LinearCodeExt<Val, EF>,
    M: Mmcs<Val>,
    L: Layout<Val, EF>,
    Challenger: FieldChallenger<Val>
        + GrindingChallenger<Witness = Val>
        + CanSampleUniformBits<Val>
        + CanSampleBits<usize>
        + CanObserve<M::Commitment>,
{
    type Val = Val;
    type Commitment = M::Commitment;
    type ProverData = BrakedownProverData<Val, EF, M, L>;
    type Proof = BrakedownProof<Val, EF, M>;
    type Error = BrakedownPcsError<M::Error>;
    type Witness = Witness<Val>;
    type OpeningProtocol = OpeningProtocol;

    fn num_vars(&self) -> usize {
        self.tensor.num_variables()
    }

    fn commit(
        &self,
        witness: Self::Witness,
        challenger: &mut Challenger,
    ) -> (Self::Commitment, Self::ProverData) {
        assert_eq!(
            witness.num_variables(),
            self.tensor.num_variables(),
            "stacked witness arity must match the tensor scheme's variable count"
        );

        // The tensor commitment and the sumcheck must operate on the *same*
        // stacked polynomial: snapshot its evaluations before consuming the
        // witness into the layout.
        let stacked = witness.poly().as_slice().to_vec();
        let (commitment, tensor) = self.tensor.commit(stacked.clone(), challenger);

        // Retain the layout purely for table access during opening.
        let layout = L::from_witness(witness);

        (
            commitment,
            BrakedownProverData {
                layout,
                stacked,
                tensor,
                _marker: PhantomData,
            },
        )
    }

    fn open(
        &self,
        _prover_data: Self::ProverData,
        _protocol: Self::OpeningProtocol,
        _challenger: &mut Challenger,
    ) -> Self::Proof {
        unimplemented!("BrakedownMlPcs supports only prescribed-point opening (open_at)")
    }

    fn verify(
        &self,
        _commitment: &Self::Commitment,
        _proof: &Self::Proof,
        _challenger: &mut Challenger,
        _protocol: Self::OpeningProtocol,
    ) -> Result<(), Self::Error> {
        unimplemented!("BrakedownMlPcs supports only prescribed-point verification (verify_at)")
    }
}

impl<EF, Val, C, M, L, Challenger> PrescribedPointPcs<EF, Challenger>
    for BrakedownMlPcs<Val, EF, C, M, L>
where
    Val: Field,
    EF: ExtensionField<Val>,
    C: LinearCodeExt<Val, EF>,
    M: Mmcs<Val>,
    L: Layout<Val, EF>,
    Challenger: FieldChallenger<Val>
        + GrindingChallenger<Witness = Val>
        + CanSampleUniformBits<Val>
        + CanSampleBits<usize>
        + CanObserve<M::Commitment>,
{
    fn open_at(
        &self,
        prover_data: Self::ProverData,
        protocol: &OpeningProtocol,
        points: &[Point<EF>],
        challenger: &mut Challenger,
    ) -> Self::Proof {
        // One prescribed point per opening batch.
        assert_eq!(protocol.num_openings(), points.len());

        let order = L::variable_order();

        // Record each batch against the shared layout verifier, computing the
        // evaluations directly from the committed tables. `add_claim_at` absorbs
        // the evals in the exact order `verify_at` replays them, so the prover
        // and verifier build the identical batched constraint.
        let mut layout_verifier = Verifier::<Val, EF>::new(&protocol.table_shapes(), L::strategy());
        let mut evals = Vec::with_capacity(points.len());
        for ((table_idx, batch), point) in protocol.iter_openings().zip(points) {
            let table = prover_data.layout.table(table_idx);
            // Direct evaluations at the prescribed (table-local) point.
            let current: Vec<EF> = batch
                .current()
                .iter()
                .map(|&col| table.poly(col).eval_base(point))
                .collect();
            // Repeat-last successor ("next row") evaluations at the same point.
            let next: Vec<EF> = batch
                .next()
                .iter()
                .map(|&col| table.poly(col).eval_next_base(point))
                .collect();
            let batch_evals = OpeningBatch::new(current, next);
            layout_verifier
                .add_claim_at(table_idx, batch, point, &batch_evals, challenger)
                .expect("prover-computed evals always match the requested batch shape");
            evals.push(batch_evals);
        }

        // Batch every recorded claim into one weight polynomial under `alpha`.
        let alpha = challenger.sample_algebra_element();
        let constraint = layout_verifier.constraint(alpha);

        // Fold all stacked variables of `P` against the batched weight `W` down
        // to a single point `r`. The product polynomial starts with zero weights;
        // `compute_sumcheck_polynomials` folds the constraint in (setting the
        // running sum to the batched claimed value) before the round loop.
        let num_variables = self.tensor.num_variables();
        let stacked_ext: Vec<EF> = prover_data.stacked.iter().map(|&v| EF::from(v)).collect();
        let product = ProductPolynomial::new_unpacked(
            order,
            Poly::new(stacked_ext),
            Poly::zero(num_variables),
        );
        let mut sumcheck_prover = SumcheckProver::new(product, EF::ZERO);
        let mut sumcheck = SumcheckData::default();
        let r = sumcheck_prover.compute_sumcheck_polynomials(
            &mut sumcheck,
            challenger,
            num_variables,
            0,
            Some(constraint),
        );

        // Tensor-open `P` at the folded point. Suffix binding evaluates the
        // stacked polynomial at the reversed challenge (matching the WHIR
        // adapter's final-value convention).
        let tensor_point = order_point(order, &r);
        let (_value, tensor) = self.tensor.open(&prover_data.tensor, &tensor_point, challenger);

        BrakedownProof {
            evals,
            sumcheck,
            tensor,
        }
    }

    fn verify_at(
        &self,
        commitment: &Self::Commitment,
        proof: &Self::Proof,
        protocol: &OpeningProtocol,
        points: &[Point<EF>],
        challenger: &mut Challenger,
    ) -> Result<Vec<OpeningEvals<EF>>, Self::Error> {
        // One prescribed point per opening batch.
        assert_eq!(protocol.num_openings(), points.len());

        if proof.evals.len() != protocol.num_openings() {
            return Err(BrakedownPcsError::OpeningCountMismatch {
                expected: protocol.num_openings(),
                got: proof.evals.len(),
            });
        }

        let order = L::variable_order();

        // Rebuild the identical batched constraint from the proof's evals.
        let mut layout_verifier = Verifier::<Val, EF>::new(&protocol.table_shapes(), L::strategy());
        for (((table_idx, batch), batch_evals), point) in
            protocol.iter_openings().zip(&proof.evals).zip(points)
        {
            layout_verifier
                .add_claim_at(table_idx, batch, point, batch_evals, challenger)
                .map_err(BrakedownPcsError::OpeningShape)?;
        }

        let alpha = challenger.sample_algebra_element();
        let constraint = layout_verifier.constraint(alpha);
        let mut claimed_eval = EF::ZERO;
        constraint.combine_evals(&mut claimed_eval);

        // Replay the sumcheck: fold every stacked variable, recovering the
        // folding point `r` and reducing `claimed_eval` to the final claim
        // `W(r) · P(r)`.
        let num_variables = self.tensor.num_variables();
        let r = p3_sumcheck::verify_final_sumcheck_rounds(
            Some(&proof.sumcheck),
            challenger,
            &mut claimed_eval,
            num_variables,
            0,
        )
        .map_err(BrakedownPcsError::Sumcheck)?;

        // Evaluate the batched weight polynomial at the folded point.
        let weight_at_r = order.eval_constraints_poly(&[constraint], &r);

        // Recover `v = P(r)` from the tensor opening (does not re-absorb the
        // commitment: the outer protocol already did).
        let tensor_point = order_point(order, &r);
        let value = self
            .tensor
            .verify_opening(commitment, &tensor_point, &proof.tensor, challenger)
            .map_err(BrakedownPcsError::Tensor)?;

        // Final consistency: the reduced claim must equal `W(r) · P(r)`.
        if claimed_eval != weight_at_r * value {
            return Err(BrakedownPcsError::FinalCheckMismatch);
        }

        Ok(proof.evals.clone())
    }
}

/// The point at which the stacked polynomial is evaluated for a given binding
/// order: prefix binding uses the folding randomness directly; suffix binding
/// reverses it (matching the WHIR adapter's final-value convention).
fn order_point<EF: Field>(order: VariableOrder, r: &Point<EF>) -> Point<EF> {
    match order {
        VariableOrder::Prefix => r.clone(),
        VariableOrder::Suffix => r.reversed(),
    }
}
