//! The `multi-stark` adapter: `QuasarMlPcs` implements [`MultilinearPcs`] and
//! [`PrescribedPointPcs`], so Quasar is a drop-in for `BrakedownMlPcs` /
//! `BasefoldMlPcs`.
//!
//! Mirrors `basefold/src/prescribed_pcs.rs`: the stacked-layout sumcheck folds
//! every table claim to a single point `r`, and the only difference is that the
//! final `P(r) = v` step is discharged by a **Quasar** opening.
//!
//! Single challenge field (`Challenge = Val`), matching this fork's ed25519
//! configuration — and Quasar's own instantiation, which needs a large field
//! for the QA distance bound anyway.

use alloc::vec::Vec;
use core::marker::PhantomData;

use p3_challenger::{
    CanObserve, CanSampleBits, CanSampleUniformBits, FieldChallenger, GrindingChallenger,
};
use p3_commit::{Mmcs, MultilinearPcs};
use p3_field::Field;
use p3_multilinear_util::point::Point;
use p3_multilinear_util::poly::Poly;
use p3_sumcheck::layout::{Layout, Table, Verifier, Witness};
use p3_sumcheck::product_polynomial::ProductPolynomial;
use p3_sumcheck::strategy::{SumcheckProver, VariableOrder};
use p3_sumcheck::{
    OpeningBatch, OpeningEvals, OpeningProtocol, PrescribedPointPcs, SumcheckData, SumcheckError,
};
use serde::{Deserialize, Serialize};

use crate::pcs::{QuasarError, QuasarPcs, QuasarProof, QuasarProverData, QuasarProverKey};

/// Quasar-backed multilinear PCS for the multi-stark opening path.
///
/// Holds the preprocessed prover key alongside the instance. Both sides can
/// build it: setup is a deterministic function of the public QA code, and the
/// verifier simply never uses the prover half.
pub struct QuasarMlPcs<Val: Field, M: Mmcs<Val>, L> {
    inner: QuasarPcs<Val, M>,
    /// Prover-only. `None` for a verifier-built instance
    /// ([`Self::new_verifier`]), which never opens anything.
    key: Option<QuasarProverKey<Val, M>>,
    _marker: PhantomData<L>,
}

impl<Val: Field, M: Mmcs<Val>, L> QuasarMlPcs<Val, M, L> {
    pub const fn new(inner: QuasarPcs<Val, M>, key: QuasarProverKey<Val, M>) -> Self {
        Self {
            inner,
            key: Some(key),
            _marker: PhantomData,
        }
    }

    /// A verifier-only instance: no preprocessed prover key, and — when
    /// `inner` came from [`QuasarPcs::new_verifier`] — no QA code either.
    ///
    /// # Panics
    ///
    /// Opening through this instance panics.
    pub const fn new_verifier(inner: QuasarPcs<Val, M>) -> Self {
        Self {
            inner,
            key: None,
            _marker: PhantomData,
        }
    }

    /// The underlying Quasar instance.
    pub const fn inner(&self) -> &QuasarPcs<Val, M> {
        &self.inner
    }
}

/// Prover handoff between commit and open.
pub struct QuasarMlProverData<Val: Field, M: Mmcs<Val>, L> {
    layout: L,
    inner: QuasarProverData<Val, M>,
    /// The stacked evaluations, kept for the layout sumcheck (the inner prover
    /// data stores them transposed for fast row combinations).
    stacked: Vec<Val>,
}

impl<Val: Field, M: Mmcs<Val>, L: Clone> Clone for QuasarMlProverData<Val, M, L>
where
    QuasarProverData<Val, M>: Clone,
{
    fn clone(&self) -> Self {
        Self {
            layout: self.layout.clone(),
            inner: self.inner.clone(),
            stacked: self.stacked.clone(),
        }
    }
}

/// Adapter opening proof: the recorded evals, the layout sumcheck, and the
/// inner Quasar `P(r) = v` proof.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(
    serialize = "Val: Serialize, M::Commitment: Serialize, M::Proof: Serialize",
    deserialize = "Val: serde::de::DeserializeOwned, M::Commitment: serde::de::DeserializeOwned, M::Proof: serde::de::DeserializeOwned"
))]
pub struct QuasarMlProof<Val: Field, M: Mmcs<Val>> {
    evals: Vec<OpeningEvals<Val>>,
    sumcheck: SumcheckData<Val, Val>,
    inner: QuasarProof<Val, M>,
}

/// Opening-verification failures.
#[derive(Debug)]
pub enum QuasarMlError<E> {
    Inner(QuasarError<E>),
    Sumcheck(SumcheckError),
    OpeningShape(SumcheckError),
    OpeningCountMismatch { expected: usize, got: usize },
    FinalCheckMismatch,
}

/// Evaluate the stacked polynomial at the folded point per binding order.
fn order_point<F: Field>(order: VariableOrder, r: &Point<F>) -> Point<F> {
    match order {
        VariableOrder::Prefix => r.clone(),
        VariableOrder::Suffix => r.reversed(),
    }
}

impl<Val, M, L, Challenger> MultilinearPcs<Val, Challenger> for QuasarMlPcs<Val, M, L>
where
    Val: Field,
    M: Mmcs<Val>,
    L: Layout<Val, Val>,
    Challenger: FieldChallenger<Val>
        + GrindingChallenger<Witness = Val>
        + CanSampleUniformBits<Val>
        + CanSampleBits<usize>
        + CanObserve<M::Commitment>
        + Clone
        + Send
        + Sync,
{
    type Val = Val;
    type Commitment = M::Commitment;
    type ProverData = QuasarMlProverData<Val, M, L>;
    type Proof = QuasarMlProof<Val, M>;
    type Error = QuasarMlError<M::Error>;
    type Witness = Witness<Val>;
    type OpeningProtocol = OpeningProtocol;

    fn num_vars(&self) -> usize {
        self.inner.num_vars()
    }

    fn commit(
        &self,
        witness: Self::Witness,
        challenger: &mut Challenger,
    ) -> (Self::Commitment, Self::ProverData) {
        assert_eq!(
            witness.num_variables(),
            self.inner.num_vars(),
            "stacked witness arity"
        );
        let stacked = witness.poly().as_slice().to_vec();
        let (commitment, inner) = self.inner.commit(stacked.clone(), challenger);
        let layout = L::from_witness(witness);
        (
            commitment,
            QuasarMlProverData {
                layout,
                inner,
                stacked,
            },
        )
    }

    fn open(
        &self,
        _pd: Self::ProverData,
        _p: Self::OpeningProtocol,
        _c: &mut Challenger,
    ) -> Self::Proof {
        unimplemented!("QuasarMlPcs supports only prescribed-point opening (open_at)")
    }

    fn verify(
        &self,
        _c: &Self::Commitment,
        _pf: &Self::Proof,
        _ch: &mut Challenger,
        _p: Self::OpeningProtocol,
    ) -> Result<(), Self::Error> {
        unimplemented!("QuasarMlPcs supports only prescribed-point verification (verify_at)")
    }
}

impl<Val, M, L, Challenger> PrescribedPointPcs<Val, Challenger> for QuasarMlPcs<Val, M, L>
where
    Val: Field,
    M: Mmcs<Val>,
    L: Layout<Val, Val>,
    Challenger: FieldChallenger<Val>
        + GrindingChallenger<Witness = Val>
        + CanSampleUniformBits<Val>
        + CanSampleBits<usize>
        + CanObserve<M::Commitment>
        + Clone
        + Send
        + Sync,
{
    fn open_at(
        &self,
        prover_data: Self::ProverData,
        protocol: &OpeningProtocol,
        points: &[Point<Val>],
        challenger: &mut Challenger,
    ) -> Self::Proof {
        assert_eq!(protocol.num_openings(), points.len());
        let order = L::variable_order();

        // Record each batch's evals at its prescribed point (absorbing them).
        let mut layout_verifier =
            Verifier::<Val, Val>::new(&protocol.table_shapes(), L::strategy());
        let mut evals = Vec::with_capacity(points.len());
        for ((table_idx, batch), point) in protocol.iter_openings().zip(points) {
            let table = prover_data.layout.table(table_idx);
            let current: Vec<Val> = batch
                .current()
                .iter()
                .map(|&col| table.poly(col).eval_base(point))
                .collect();
            let next: Vec<Val> = batch
                .next()
                .iter()
                .map(|&col| table.poly(col).eval_next_base(point))
                .collect();
            let batch_evals = OpeningBatch::new(current, next);
            layout_verifier
                .add_claim_at(table_idx, batch, point, &batch_evals, challenger)
                .expect("prover-computed evals match the requested shape");
            evals.push(batch_evals);
        }

        // Batch under alpha, fold all stacked variables to a point r.
        let alpha = challenger.sample_algebra_element();
        let constraint = layout_verifier.constraint(alpha);
        let num_variables = self.inner.num_vars();
        let product = ProductPolynomial::new_unpacked(
            order,
            Poly::new(prover_data.stacked.clone()),
            Poly::zero(num_variables),
        );
        let mut sumcheck_prover = SumcheckProver::new(product, Val::ZERO);
        let mut sumcheck = SumcheckData::default();
        let r = sumcheck_prover.compute_sumcheck_polynomials(
            &mut sumcheck,
            challenger,
            num_variables,
            0,
            Some(constraint),
        );

        // Discharge P(r) = v with a Quasar opening at the folded point.
        let z = order_point(order, &r);
        let (_v, inner) = self
            .inner
            .open(
                self.key
                    .as_ref()
                    .expect("this is a verifier-only QuasarMlPcs: it cannot open"),
                &prover_data.inner,
                z.as_slice(),
                challenger,
            );

        QuasarMlProof {
            evals,
            sumcheck,
            inner,
        }
    }

    fn verify_at(
        &self,
        commitment: &Self::Commitment,
        proof: &Self::Proof,
        protocol: &OpeningProtocol,
        points: &[Point<Val>],
        challenger: &mut Challenger,
    ) -> Result<Vec<OpeningEvals<Val>>, Self::Error> {
        assert_eq!(protocol.num_openings(), points.len());
        if proof.evals.len() != protocol.num_openings() {
            return Err(QuasarMlError::OpeningCountMismatch {
                expected: protocol.num_openings(),
                got: proof.evals.len(),
            });
        }
        let order = L::variable_order();

        let mut layout_verifier =
            Verifier::<Val, Val>::new(&protocol.table_shapes(), L::strategy());
        for (((table_idx, batch), batch_evals), point) in
            protocol.iter_openings().zip(&proof.evals).zip(points)
        {
            layout_verifier
                .add_claim_at(table_idx, batch, point, batch_evals, challenger)
                .map_err(QuasarMlError::OpeningShape)?;
        }

        let alpha = challenger.sample_algebra_element();
        let constraint = layout_verifier.constraint(alpha);
        let mut claimed_eval = Val::ZERO;
        constraint.combine_evals(&mut claimed_eval);

        let num_variables = self.inner.num_vars();
        let r = p3_sumcheck::verify_final_sumcheck_rounds(
            Some(&proof.sumcheck),
            challenger,
            &mut claimed_eval,
            num_variables,
            0,
        )
        .map_err(QuasarMlError::Sumcheck)?;

        let weight_at_r = order.eval_constraints_poly(&[constraint], &r);

        let z = order_point(order, &r);
        let value = self
            .inner
            .verify_opening(commitment, z.as_slice(), &proof.inner, challenger)
            .map_err(QuasarMlError::Inner)?;

        if claimed_eval != weight_at_r * value {
            return Err(QuasarMlError::FinalCheckMismatch);
        }
        Ok(proof.evals.clone())
    }
}

impl<Val: Field, M: Mmcs<Val>, L: Layout<Val, Val>> QuasarMlProverData<Val, M, L> {
    /// Table access for `MultiStarkConfig::committed_table`.
    pub fn table(&self, index: usize) -> &Table<Val> {
        self.layout.table(index)
    }
}
