//! M3: the multi-stark adapter — `BasefoldMlPcs` implements `MultilinearPcs` + `PrescribedPointPcs`
//! so Basefold is a drop-in for `BrakedownMlPcs`. Mirrors `brakedown/src/prescribed_pcs.rs`: the
//! stacked-layout sumcheck folds all table claims to a point `r`; the only difference is the final
//! `P(r)=v` step is discharged by the **Basefold** opening (`merkle_pcs`) instead of the tensor one.
//!
//! Single challenge field (`Challenge = Val`), matching alp's eddsa config (`EF = Ed25519 = F`).

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

use crate::merkle_pcs::{BasefoldPcs, BasefoldProof, BasefoldProverData};

/// Basefold-backed multilinear PCS for the multi-stark opening path.
pub struct BasefoldMlPcs<Val: Field, M: Mmcs<Val>, L> {
    inner: BasefoldPcs<Val, M>,
    _marker: PhantomData<L>,
}

impl<Val: Field, M: Mmcs<Val>, L> BasefoldMlPcs<Val, M, L> {
    pub const fn new(inner: BasefoldPcs<Val, M>) -> Self {
        Self { inner, _marker: PhantomData }
    }
}

/// Prover handoff between commit and open.
pub struct BasefoldMlProverData<Val: Field, M: Mmcs<Val>, L> {
    layout: L,
    inner: BasefoldProverData<Val, M>,
}

impl<Val: Field, M: Mmcs<Val>, L: Clone> Clone for BasefoldMlProverData<Val, M, L>
where
    BasefoldProverData<Val, M>: Clone,
{
    fn clone(&self) -> Self {
        Self { layout: self.layout.clone(), inner: self.inner.clone() }
    }
}

/// Adapter opening proof: recorded evals + the layout sumcheck + the inner Basefold `P(r)=v` proof.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(
    serialize = "Val: Serialize, M::Commitment: Serialize, M::Proof: Serialize",
    deserialize = "Val: serde::de::DeserializeOwned, M::Commitment: serde::de::DeserializeOwned, M::Proof: serde::de::DeserializeOwned"
))]
pub struct BasefoldMlProof<Val: Field, M: Mmcs<Val>> {
    evals: Vec<OpeningEvals<Val>>,
    sumcheck: SumcheckData<Val, Val>,
    inner: BasefoldProof<Val, M>,
}

/// Opening-verification failures.
#[derive(Debug)]
pub enum BasefoldMlError {
    Inner(&'static str),
    Sumcheck(SumcheckError),
    OpeningShape(SumcheckError),
    OpeningCountMismatch { expected: usize, got: usize },
    FinalCheckMismatch,
}

/// Evaluate the stacked polynomial at the folded point per binding order (mirrors brakedown).
fn order_point<F: Field>(order: VariableOrder, r: &Point<F>) -> Point<F> {
    match order {
        VariableOrder::Prefix => r.clone(),
        VariableOrder::Suffix => r.reversed(),
    }
}

impl<Val, M, L, Challenger> MultilinearPcs<Val, Challenger> for BasefoldMlPcs<Val, M, L>
where
    Val: Field,
    M: Mmcs<Val>,
    L: Layout<Val, Val>,
    Challenger: FieldChallenger<Val>
        + GrindingChallenger<Witness = Val>
        + CanSampleUniformBits<Val>
        + CanSampleBits<usize>
        + CanObserve<M::Commitment>,
{
    type Val = Val;
    type Commitment = M::Commitment;
    type ProverData = BasefoldMlProverData<Val, M, L>;
    type Proof = BasefoldMlProof<Val, M>;
    type Error = BasefoldMlError;
    type Witness = Witness<Val>;
    type OpeningProtocol = OpeningProtocol;

    fn num_vars(&self) -> usize {
        self.inner.num_vars()
    }

    fn commit(&self, witness: Self::Witness, challenger: &mut Challenger) -> (Self::Commitment, Self::ProverData) {
        assert_eq!(witness.num_variables(), self.inner.num_vars(), "stacked witness arity");
        let stacked = witness.poly().as_slice().to_vec();
        let (commitment, inner) = self.inner.commit(stacked, challenger);
        let layout = L::from_witness(witness);
        (commitment, BasefoldMlProverData { layout, inner })
    }

    fn open(&self, _pd: Self::ProverData, _p: Self::OpeningProtocol, _c: &mut Challenger) -> Self::Proof {
        unimplemented!("BasefoldMlPcs supports only prescribed-point opening (open_at)")
    }
    fn verify(&self, _c: &Self::Commitment, _pf: &Self::Proof, _ch: &mut Challenger, _p: Self::OpeningProtocol) -> Result<(), Self::Error> {
        unimplemented!("BasefoldMlPcs supports only prescribed-point verification (verify_at)")
    }
}

impl<Val, M, L, Challenger> PrescribedPointPcs<Val, Challenger> for BasefoldMlPcs<Val, M, L>
where
    Val: Field,
    M: Mmcs<Val>,
    L: Layout<Val, Val>,
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
        points: &[Point<Val>],
        challenger: &mut Challenger,
    ) -> Self::Proof {
        assert_eq!(protocol.num_openings(), points.len());
        let order = L::variable_order();

        // Record each batch's evals at its prescribed point (absorbs them into the transcript).
        let mut layout_verifier = Verifier::<Val, Val>::new(&protocol.table_shapes(), L::strategy());
        let mut evals = Vec::with_capacity(points.len());
        for ((table_idx, batch), point) in protocol.iter_openings().zip(points) {
            let table = prover_data.layout.table(table_idx);
            let current: Vec<Val> = batch.current().iter().map(|&col| table.poly(col).eval_base(point)).collect();
            let next: Vec<Val> = batch.next().iter().map(|&col| table.poly(col).eval_next_base(point)).collect();
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
            Poly::new(prover_data.inner_evals().to_vec()),
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

        // Discharge P(r) = v via the Basefold opening at the folded point.
        let z = order_point(order, &r);
        let (_v, inner) = self.inner.open(&prover_data.inner, z.as_slice(), challenger);

        BasefoldMlProof { evals, sumcheck, inner }
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
            return Err(BasefoldMlError::OpeningCountMismatch {
                expected: protocol.num_openings(),
                got: proof.evals.len(),
            });
        }
        let order = L::variable_order();

        let mut layout_verifier = Verifier::<Val, Val>::new(&protocol.table_shapes(), L::strategy());
        for (((table_idx, batch), batch_evals), point) in protocol.iter_openings().zip(&proof.evals).zip(points) {
            layout_verifier
                .add_claim_at(table_idx, batch, point, batch_evals, challenger)
                .map_err(BasefoldMlError::OpeningShape)?;
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
        .map_err(BasefoldMlError::Sumcheck)?;

        let weight_at_r = order.eval_constraints_poly(&[constraint], &r);

        let z = order_point(order, &r);
        let value = self
            .inner
            .verify_opening(commitment, z.as_slice(), &proof.inner, challenger)
            .map_err(BasefoldMlError::Inner)?;

        if claimed_eval != weight_at_r * value {
            return Err(BasefoldMlError::FinalCheckMismatch);
        }
        Ok(proof.evals.clone())
    }
}

impl<Val: Field, M: Mmcs<Val>, L> BasefoldMlProverData<Val, M, L> {
    fn inner_evals(&self) -> &[Val] {
        self.inner.evals_slice()
    }
}

impl<Val: Field, M: Mmcs<Val>, L: Layout<Val, Val>> BasefoldMlProverData<Val, M, L> {
    /// Table access for `MultiStarkConfig::committed_table`.
    pub fn table(&self, index: usize) -> &Table<Val> {
        self.layout.table(index)
    }
}
