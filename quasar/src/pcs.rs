//! The Quasar polynomial commitment scheme.

use alloc::vec;
use alloc::vec::Vec;

use p3_basefold::{BasefoldBatchProverData, BasefoldPcs, BasefoldProof};
use p3_challenger::{CanObserve, CanSampleBits, FieldChallenger, GrindingChallenger};
use p3_commit::Mmcs;
use p3_field::Field;
use p3_matrix::Dimensions;
use p3_matrix::dense::RowMajorMatrix;
use p3_maybe_rayon::prelude::*;
use p3_multilinear_util::poly::Poly;
use p3_brakedown::LinearCode;
use p3_qa_code::QaCode;
use p3_sumcheck::generic_degree::{GenericDegreeProof, RoundProver};
use p3_util::log2_ceil_usize;
use serde::{Deserialize, Serialize};

use crate::params::QuasarParams;
use crate::sumchecks::{Deg2Prover, Deg3Prover};
use crate::util::{dot, eq_at_index, eq_eval, eq_table, eval_mle, powers, wht_mle_eval, wht_mle_table};

/// Which committed vector an evaluation claim refers to.
///
/// Online vectors are indexed `0..2c`: `0..c` are the blocks of the proximity
/// codeword `u` (so `0` is the message block `v`), and `c..2c` those of the
/// evaluation codeword `q` (so `c` is `p`). Coefficient vectors are the `c − 1`
/// preprocessed `2^{-n₂}a_i'`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Target {
    Online(usize),
    Coeff(usize),
}

/// One evaluation claim to be discharged by a batched Basefold opening.
struct Claim<F> {
    target: Target,
    point: Vec<F>,
    value: F,
}

/// Verification failures.
#[derive(Debug)]
pub enum QuasarError<E> {
    /// A committed-column Merkle path did not verify.
    Merkle(E),
    /// A Basefold opening failed.
    Basefold(&'static str),
    /// A sumcheck transcript was malformed.
    Sumcheck(&'static str),
    /// The proof carried the wrong number of openings or claimed values.
    MalformedProof(&'static str),
    /// The selector sumcheck's claimed sum disagreed with the opened columns.
    ColumnConsistency,
    /// The selector sumcheck's terminal claim disagreed with the block
    /// evaluations.
    SelectorTerminal,
    /// The QA-encoding claim disagreed with the parity-block evaluations.
    EncodingClaim,
    /// The outer encoding sumcheck's terminal claim failed.
    EncodingOuterTerminal,
    /// The inner encoding sumcheck's terminal claim failed.
    EncodingInnerTerminal,
    /// The claim-reduction sumcheck disagreed with the collected claims.
    ReductionClaim,
    /// The claim-reduction sumcheck's terminal claim failed.
    ReductionTerminal,
    /// A batched Basefold opening bound a different value than the reduction
    /// left standing.
    OpeningValue,
    /// The proof's claimed evaluation is not the one the caller asked about.
    ClaimedValueMismatch,
    /// The proof-of-work witness gating the column sample did not check out.
    InvalidProofOfWork,
}

/// A Quasar evaluation proof.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(
    serialize = "F: Serialize, M::Commitment: Serialize, M::Proof: Serialize",
    deserialize = "F: serde::de::DeserializeOwned, M::Commitment: serde::de::DeserializeOwned, M::Proof: serde::de::DeserializeOwned"
))]
pub struct QuasarProof<F: Field, M: Mmcs<F>> {
    /// The evaluation this proof binds the commitment to. `verify` checks it
    /// against the caller's claim; `verify_opening` hands it to an outer
    /// protocol to check against its own.
    pub claimed_value: F,
    /// Batched Basefold commitment to the `2c` online blocks.
    pub blocks_commitment: M::Commitment,
    /// Proof-of-work witness bound just before the columns were sampled.
    /// `F::ZERO` when grinding is off.
    pub pow_witness: F,
    /// The sampled columns of the committed codeword matrix, under ONE pruned
    /// authentication path: `column_openings[k][0]` is column `indices[k]`.
    pub column_openings: Vec<Vec<Vec<F>>>,
    /// The pruned paths covering every sampled column.
    pub column_proof: M::MultiProof,
    /// Sampled-column consistency sumcheck, over `log(c·N₂)` variables.
    pub selector: GenericDegreeProof<F, F>,
    /// The `2c` online blocks at the selector's terminal point.
    pub block_evals_at_e1: Vec<F>,
    /// The `2(c−1)` parity blocks at `γ`.
    pub parity_evals_at_gamma: Vec<F>,
    /// Outer QA-encoding sumcheck (degree 3, `n₂` variables).
    pub encoding_outer: GenericDegreeProof<F, F>,
    /// The `c−1` preprocessed coefficient vectors at the outer terminal point.
    pub coeff_evals_at_e2: Vec<F>,
    /// Inner QA-encoding sumcheck (degree 2, `n₂` variables).
    pub encoding_inner: GenericDegreeProof<F, F>,
    /// `[v(e₃), p(e₃)]` — the two committed input blocks at the inner terminal
    /// point.
    pub input_evals_at_e3: [F; 2],
    /// Claim-reduction sumcheck collapsing every claim to one point `ζ`.
    pub reduction: GenericDegreeProof<F, F>,
    /// Every committed vector at `ζ`: the `2c` online blocks then the `c−1`
    /// coefficient vectors.
    pub evals_at_zeta: Vec<F>,
    /// One batched Basefold opening of the online blocks at `ζ`.
    pub online_opening: BasefoldProof<F, M>,
    /// One batched Basefold opening of the preprocessed coefficients at `ζ`.
    pub coeff_opening: BasefoldProof<F, M>,
}

/// Prover state retained between commit and open.
pub struct QuasarProverData<F: Field, M: Mmcs<F>> {
    /// The message matrix in `N₂ × N₁` layout (`msg[j][i] = M[i][j]`), so a
    /// row-combination is a contiguous dot product.
    msg: RowMajorMatrix<F>,
    mmcs_data: M::ProverData<RowMajorMatrix<F>>,
}

// The multi-stark prover clones prover data across its opening rounds, so this
// is required by the `PrescribedPointPcs` path (mirrors `BasefoldProverData`).
impl<F: Field, M: Mmcs<F>> Clone for QuasarProverData<F, M>
where
    M::ProverData<RowMajorMatrix<F>>: Clone,
{
    fn clone(&self) -> Self {
        Self {
            msg: self.msg.clone(),
            mmcs_data: self.mmcs_data.clone(),
        }
    }
}

/// Preprocessed prover material: the committed QA coefficient vectors.
pub struct QuasarProverKey<F: Field, M: Mmcs<F>> {
    coeff_data: BasefoldBatchProverData<F, M>,
}

/// Quasar over a single large prime field.
///
/// Challenges are drawn from `F` itself rather than an extension — as in the
/// paper (Mersenne127) and this fork's ed25519 configuration, where the field is
/// already wide enough for the Fiat–Shamir terms. QA codes need a large field
/// for their distance bound anyway, so this is not a restriction in practice.
pub struct QuasarPcs<F: Field, M: Mmcs<F>> {
    params: QuasarParams,
    /// The QA code. Prover-only: every use is in [`Self::commit`] and
    /// [`Self::open`], so a verifier-only instance ([`Self::new_verifier`])
    /// leaves this `None` rather than sampling `(c−1)·N₂` field elements it
    /// would never read.
    code: Option<QaCode<F>>,
    mmcs: M,
    basefold: BasefoldPcs<F, M>,
    coeff_commitment: M::Commitment,
    /// Proof-of-work bits ground immediately before the committed columns are
    /// sampled. Each bit buys one bit of column-sampling soundness, so
    /// `params.num_queries` can be derived from `λ − pow_bits` instead of `λ`.
    ///
    /// The inner Basefold instance carries its own `pow_bits` for its query
    /// phase; the two are independent.
    pow_bits: usize,
}

impl<F: Field, M: Mmcs<F>> QuasarPcs<F, M> {
    /// Set up an instance, committing the public QA coefficient vectors.
    ///
    /// `basefold` must be configured for length-`N₂` messages: those are the
    /// vectors code switching moves the claims onto.
    ///
    /// # Returns
    ///
    /// The instance (usable by prover and verifier) and the prover-only
    /// preprocessed data.
    pub fn new(
        params: QuasarParams,
        code: QaCode<F>,
        mmcs: M,
        basefold: BasefoldPcs<F, M>,
    ) -> (Self, QuasarProverKey<F, M>) {
        assert_eq!(
            code.log_message_len(),
            params.log_cols,
            "QA code message length must be N2"
        );
        assert_eq!(code.index(), params.index, "QA index mismatch");
        assert_eq!(
            basefold.num_vars(),
            params.log_cols,
            "Basefold must commit length-N2 vectors"
        );
        assert!(
            params.num_queries <= params.codeword_len(),
            "cannot open more columns than the codeword has"
        );

        let (coeff_commitment, coeff_data) =
            basefold.commit_many_silent(code.scaled_coeffs().to_vec());

        (
            Self {
                params,
                code: Some(code),
                mmcs,
                basefold,
                coeff_commitment,
                pow_bits: 0,
            },
            QuasarProverKey { coeff_data },
        )
    }

    /// A **verifier-only** instance, given the coefficient commitment instead
    /// of the QA code that produces it.
    ///
    /// [`Self::new`] ends with `basefold.commit_many_silent(code.scaled_coeffs())`
    /// — a full encode-and-Merkle-commit of the `c−1` public coefficient
    /// vectors. The verifier consumes only the resulting digest, and the
    /// digest is a *public parameter*: it is fixed by the code seed, exactly
    /// like the preprocessed commitment in a verifying key. Recomputing it
    /// costs `O(N₂)` field work and hashing, which is free on a host that sets
    /// up once and expensive in a recursion circuit that sets up per proof.
    ///
    /// Pass the commitment `new` returned (`Self::coeff_commitment`).
    ///
    /// # Panics
    ///
    /// Calling [`Self::commit`] or [`Self::open`] on the result panics: there
    /// is no QA code to encode with.
    pub fn new_verifier(
        params: QuasarParams,
        mmcs: M,
        basefold: BasefoldPcs<F, M>,
        coeff_commitment: M::Commitment,
    ) -> Self {
        assert_eq!(
            basefold.num_vars(),
            params.log_cols,
            "Basefold must commit length-N2 vectors"
        );
        assert!(
            params.num_queries <= params.codeword_len(),
            "cannot open more columns than the codeword has"
        );
        Self { params, code: None, mmcs, basefold, coeff_commitment, pow_bits: 0 }
    }

    /// The coefficient commitment — a public parameter, to be handed to
    /// [`Self::new_verifier`].
    pub const fn coeff_commitment(&self) -> &M::Commitment {
        &self.coeff_commitment
    }

    /// The QA code, for the prover paths that encode with it.
    ///
    /// # Panics
    ///
    /// If this is a verifier-only instance.
    fn code(&self) -> &QaCode<F> {
        self.code
            .as_ref()
            .expect("this is a verifier-only QuasarPcs: it cannot commit or open")
    }

    /// Add `pow_bits` of proof of work before the column-sampling phase.
    ///
    /// The caller is responsible for having derived `params.num_queries` from
    /// the *reduced* target `λ − pow_bits`; this only installs the gate.
    #[must_use]
    pub const fn with_grinding(mut self, pow_bits: usize) -> Self {
        self.pow_bits = pow_bits;
        self
    }

    /// The instance's shape.
    pub const fn params(&self) -> &QuasarParams {
        &self.params
    }

    /// Number of variables of the committed polynomial.
    pub const fn num_vars(&self) -> usize {
        self.params.num_vars()
    }

    /// Commit to a multilinear polynomial given by its `2^m` evaluations.
    ///
    /// Reshape to `N₁ × N₂`, QA-encode every row, Merkle-commit the columns of
    /// the encoded matrix. Identical in shape to a Brakedown-style tensor
    /// commitment; only the row code differs.
    pub fn commit<Ch>(
        &self,
        evals: Vec<F>,
        challenger: &mut Ch,
    ) -> (M::Commitment, QuasarProverData<F, M>)
    where
        Ch: CanObserve<M::Commitment>,
    {
        let (rows, cols) = (self.params.num_rows(), self.params.num_cols());
        assert_eq!(evals.len(), rows * cols, "eval count must be N1·N2");

        // Row-major M (rows × cols), transposed so the code — which encodes
        // along its height — encodes each row of M.
        let msg = tracing::info_span!("quasar_reshape")
            .in_scope(|| RowMajorMatrix::new(evals, cols).transpose()); // cols × rows
        let encoded =
            tracing::info_span!("quasar_qa_encode").in_scope(|| self.code().encode(&msg)); // Ĉᵗ
        let (commitment, mmcs_data) =
            tracing::info_span!("quasar_column_merkle").in_scope(|| self.mmcs.commit_matrix(encoded));
        challenger.observe(commitment.clone());
        (commitment, QuasarProverData { msg, mmcs_data })
    }

    /// Prove `f(point) = y`, returning `y` and the proof.
    pub fn open<Ch>(
        &self,
        pk: &QuasarProverKey<F, M>,
        pd: &QuasarProverData<F, M>,
        point: &[F],
        challenger: &mut Ch,
    ) -> (F, QuasarProof<F, M>)
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>
            + Clone
            + Send
            + Sync,
    {
        let params = &self.params;
        assert_eq!(point.len(), params.num_vars(), "point arity");
        let (n2, c) = (params.log_cols, params.index);
        let cols = params.num_cols();

        challenger.observe(self.coeff_commitment.clone());

        // --- Row combinations: the proximity row and the evaluation row. ---
        let r: Vec<F> = (0..params.num_rows())
            .map(|_| challenger.sample_algebra_element())
            .collect();
        let (z_hi, z_lo) = point.split_at(params.log_rows);
        let z_left = eq_table(z_hi);
        let z_right = eq_table(z_lo);

        let (v, p_row) = tracing::info_span!("quasar_row_combine").in_scope(|| {
            (
                self.combine_rows(&pd.msg, &r),
                self.combine_rows(&pd.msg, &z_left),
            )
        });
        let value = dot(&p_row, &z_right);

        // --- Code switching: encode both rows, commit their blocks. ---
        let (u, q) = tracing::info_span!("quasar_row_encode")
            .in_scope(|| (self.code().encode_vec(&v), self.code().encode_vec(&p_row)));
        let blocks: Vec<Vec<F>> = (0..c)
            .map(|i| u[i * cols..(i + 1) * cols].to_vec())
            .chain((0..c).map(|i| q[i * cols..(i + 1) * cols].to_vec()))
            .collect();
        let (blocks_commitment, blocks_data) = tracing::info_span!("quasar_block_commit")
            .in_scope(|| self.basefold.commit_many(blocks, challenger));

        // --- Sampled columns of the committed codeword matrix. ---
        // Grind first, so a prover who dislikes its column set must redo the
        // work to draw another.
        let pow_witness = tracing::info_span!("quasar_grind")
            .in_scope(|| p3_basefold::grind_parallel(challenger, self.pow_bits));
        let indices = self.sample_columns(challenger);
        let (column_openings, column_proof) = tracing::info_span!("quasar_column_open")
            .in_scope(|| self.mmcs.open_multi_batch(&indices, &pd.mmcs_data));

        // --- Selector sumcheck: the opened columns agree with u and q. ---
        let tau: F = challenger.sample_algebra_element();
        let sigma: F = challenger.sample_algebra_element();
        let selector_table = self.selector_table(&indices, tau);
        let combined: Vec<F> = u.par_iter().zip(q.par_iter()).map(|(&a, &b)| a + sigma * b).collect();

        let mut selector_prover =
            Deg2Prover::new(vec![(Poly::new(selector_table), Poly::new(combined))]);
        let selector_claim = selector_prover.claimed_sum();
        challenger.observe_algebra_element(selector_claim);
        let (selector, e1) = tracing::info_span!("quasar_selector_sumcheck").in_scope(|| {
            selector_prover.prove::<F, _>(challenger, params.log_codeword_len(), 2, 0, selector_claim)
        });

        let e1_lo = e1.as_slice()[params.log_index()..].to_vec();
        let block_evals_at_e1: Vec<F> = (0..params.num_online_blocks())
            .map(|j| eval_mle(blocks_data.poly(j), &e1_lo))
            .collect();
        challenger.observe_algebra_slice(&block_evals_at_e1);

        // --- QA-encoding check, outer layer. ---
        let gamma: Vec<F> = (0..n2).map(|_| challenger.sample_algebra_element()).collect();
        let alpha: F = challenger.sample_algebra_element();
        let beta: F = challenger.sample_algebra_element();

        let parity_evals_at_gamma: Vec<F> = (1..c)
            .map(|i| eval_mle(blocks_data.poly(i), &gamma))
            .chain((1..c).map(|i| eval_mle(blocks_data.poly(c + i), &gamma)))
            .collect();
        challenger.observe_algebra_slice(&parity_evals_at_gamma);

        // The three factors: H̃(γ,·), the α-batched public coefficients, and the
        // β-batched encoding witness WHT(v) + β·WHT(p) — which is never
        // committed (that is optimization 2 of the paper).
        let coeff_combined = self.combine_coeffs(alpha);
        let mid: Vec<F> = {
            let vm = self.code().middle_vec(&v);
            let pm = self.code().middle_vec(&p_row);
            vm.par_iter().zip(pm.par_iter()).map(|(&a, &b)| a + beta * b).collect()
        };
        let mut outer_prover = Deg3Prover::new(vec![(
            Poly::new(wht_mle_table(&gamma)),
            Poly::new(coeff_combined),
            Poly::new(mid),
        )]);
        let outer_claim = outer_prover.claimed_sum();
        debug_assert_eq!(
            outer_claim,
            batched_parity_claim(c, alpha, beta, &parity_evals_at_gamma),
            "outer sumcheck claim must be the batched parity-block evaluations"
        );
        challenger.observe_algebra_element(outer_claim);
        let (encoding_outer, e2) = tracing::info_span!("quasar_encoding_outer")
            .in_scope(|| outer_prover.prove::<F, _>(challenger, n2, 3, 0, outer_claim));
        let e2 = e2.as_slice().to_vec();

        let coeff_evals_at_e2: Vec<F> = self
            .code()
            .scaled_coeffs()
            .iter()
            .map(|e| eval_mle(e, &e2))
            .collect();
        challenger.observe_algebra_slice(&coeff_evals_at_e2);

        // --- QA-encoding check, inner layer: reduce WHT(v)(e2) to v(e3). ---
        let input: Vec<F> = v.par_iter().zip(p_row.par_iter()).map(|(&a, &b)| a + beta * b).collect();
        let mut inner_prover =
            Deg2Prover::new(vec![(Poly::new(wht_mle_table(&e2)), Poly::new(input))]);
        let inner_claim = inner_prover.claimed_sum();
        challenger.observe_algebra_element(inner_claim);
        let (encoding_inner, e3) = tracing::info_span!("quasar_encoding_inner")
            .in_scope(|| inner_prover.prove::<F, _>(challenger, n2, 2, 0, inner_claim));
        let e3 = e3.as_slice().to_vec();

        let input_evals_at_e3 = [eval_mle(&v, &e3), eval_mle(&p_row, &e3)];
        challenger.observe_algebra_slice(&input_evals_at_e3);

        // --- Collapse every claim to one point, then open twice. ---
        let claims = self.collect_claims(
            &e1_lo,
            &gamma,
            &e2,
            &e3,
            z_lo,
            value,
            &block_evals_at_e1,
            &parity_evals_at_gamma,
            &input_evals_at_e3,
            &coeff_evals_at_e2,
        );

        let rho: F = challenger.sample_algebra_element();
        let weights = self.claim_weights(&claims, rho);
        let pairs: Vec<(Poly<F>, Poly<F>)> = weights
            .online
            .into_iter()
            .enumerate()
            .map(|(j, w)| (Poly::new(w), Poly::new(blocks_data.poly(j).to_vec())))
            .chain(
                weights
                    .coeff
                    .into_iter()
                    .enumerate()
                    .map(|(i, w)| (Poly::new(w), Poly::new(self.code().scaled_coeffs()[i].clone()))),
            )
            .collect();

        let mut reduction_prover = Deg2Prover::new(pairs);
        let reduction_claim = reduction_prover.claimed_sum();
        debug_assert_eq!(
            reduction_claim,
            batched_claim_value(&claims, rho),
            "reduction claim must be the rho-combination of the claims"
        );
        challenger.observe_algebra_element(reduction_claim);
        let (reduction, zeta) = tracing::info_span!("quasar_claim_reduction")
            .in_scope(|| reduction_prover.prove::<F, _>(challenger, n2, 2, 0, reduction_claim));
        let zeta = zeta.as_slice().to_vec();

        let evals_at_zeta: Vec<F> = (0..params.num_online_blocks())
            .map(|j| eval_mle(blocks_data.poly(j), &zeta))
            .chain(
                self.code()
                    .scaled_coeffs()
                    .iter()
                    .map(|e| eval_mle(e, &zeta)),
            )
            .collect();
        challenger.observe_algebra_slice(&evals_at_zeta);

        let eta: F = challenger.sample_algebra_element();
        let (_, online_opening) = tracing::info_span!("quasar_basefold_open_online").in_scope(|| {
            self.basefold.open_batch(
                &blocks_data,
                &powers(eta, params.num_online_blocks()),
                &zeta,
                challenger,
            )
        });
        let (_, coeff_opening) = tracing::info_span!("quasar_basefold_open_coeff").in_scope(|| {
            self.basefold.open_batch(
                &pk.coeff_data,
                &powers(eta, params.num_coeff_vectors()),
                &zeta,
                challenger,
            )
        });

        (
            value,
            QuasarProof {
                claimed_value: value,
                blocks_commitment,
                pow_witness,
                column_openings,
                column_proof,
                selector,
                block_evals_at_e1,
                parity_evals_at_gamma,
                encoding_outer,
                coeff_evals_at_e2,
                encoding_inner,
                input_evals_at_e3,
                reduction,
                evals_at_zeta,
                online_opening,
                coeff_opening,
            },
        )
    }

    /// Verify that the commitment opens to `value` at `point`.
    pub fn verify<Ch>(
        &self,
        commitment: &M::Commitment,
        point: &[F],
        value: F,
        proof: &QuasarProof<F, M>,
        challenger: &mut Ch,
    ) -> Result<(), QuasarError<M::Error>>
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>
            + Clone
            + Send
            + Sync,
    {
        if proof.claimed_value != value {
            return Err(QuasarError::ClaimedValueMismatch);
        }
        self.verify_core(true, commitment, point, proof, challenger)
    }

    /// Adapter-facing verify: recover the value this opening binds `point` to,
    /// for an outer protocol to check against its own batched claim, and do
    /// **not** re-absorb the commitment (the outer verifier already did).
    ///
    /// Mirrors `BasefoldPcs::verify_opening` and `TensorPcs::verify_opening`.
    pub fn verify_opening<Ch>(
        &self,
        commitment: &M::Commitment,
        point: &[F],
        proof: &QuasarProof<F, M>,
        challenger: &mut Ch,
    ) -> Result<F, QuasarError<M::Error>>
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>
            + Clone
            + Send
            + Sync,
    {
        self.verify_core(false, commitment, point, proof, challenger)?;
        Ok(proof.claimed_value)
    }

    /// The verification engine. `observe_commitment = false` for the adapter
    /// path.
    fn verify_core<Ch>(
        &self,
        observe_commitment: bool,
        commitment: &M::Commitment,
        point: &[F],
        proof: &QuasarProof<F, M>,
        challenger: &mut Ch,
    ) -> Result<(), QuasarError<M::Error>>
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>
            + Clone
            + Send
            + Sync,
    {
        let value = proof.claimed_value;
        let params = &self.params;
        assert_eq!(point.len(), params.num_vars(), "point arity");
        let (n2, c) = (params.log_cols, params.index);

        if proof.block_evals_at_e1.len() != params.num_online_blocks()
            || proof.parity_evals_at_gamma.len() != 2 * (c - 1)
            || proof.coeff_evals_at_e2.len() != params.num_coeff_vectors()
            || proof.evals_at_zeta.len() != params.num_online_blocks() + params.num_coeff_vectors()
        {
            return Err(QuasarError::MalformedProof("claimed-value count"));
        }

        // Replay the transcript exactly as the prover drove it.
        if observe_commitment {
            challenger.observe(commitment.clone());
        }
        challenger.observe(self.coeff_commitment.clone());

        let r: Vec<F> = (0..params.num_rows())
            .map(|_| challenger.sample_algebra_element())
            .collect();
        let (z_hi, z_lo) = point.split_at(params.log_rows);
        let z_left = eq_table(z_hi);

        challenger.observe(proof.blocks_commitment.clone());
        if self.pow_bits > 0 && !challenger.check_witness(self.pow_bits, proof.pow_witness) {
            return Err(QuasarError::InvalidProofOfWork);
        }
        let indices = self.sample_columns(challenger);
        if proof.column_openings.len() != indices.len() {
            return Err(QuasarError::MalformedProof("column opening count"));
        }

        let tau: F = challenger.sample_algebra_element();
        let sigma: F = challenger.sample_algebra_element();

        // Check the column paths and rebuild the selector claim from them.
        // One fold with `r + σ·z_left` covers both branches at once.
        let fold_weights: Vec<F> = r
            .iter()
            .zip(&z_left)
            .map(|(&a, &b)| a + sigma * b)
            .collect();
        let dims = [Dimensions {
            width: params.num_rows(),
            height: params.codeword_len(),
        }];
        if proof
            .column_openings
            .iter()
            .any(|v| v.len() != 1 || v[0].len() != params.num_rows())
        {
            return Err(QuasarError::MalformedProof("column height"));
        }
        self.mmcs
            .verify_multi_batch(
                commitment,
                &dims,
                &indices,
                &proof.column_openings,
                &proof.column_proof,
            )
            .map_err(QuasarError::Merkle)?;

        let tau_powers = powers(tau, indices.len());
        let mut expected_selector_claim = F::ZERO;
        for (k, opening) in proof.column_openings.iter().enumerate() {
            expected_selector_claim += tau_powers[k] * dot(&fold_weights, &opening[0]);
        }
        if proof.selector.claimed_sum != expected_selector_claim {
            return Err(QuasarError::ColumnConsistency);
        }

        // COST ABLATION ONLY — never enable in a real verifier. Stops after the
        // sampled-column Merkle check and the selector claim, so a recursion
        // guest's instruction count isolates the column phase.
        #[cfg(feature = "ablate-after-columns")]
        return Ok(());

        #[allow(unreachable_code)]
        let (e1, selector_final) = proof
            .selector
            .verify(challenger, params.log_codeword_len(), 2, 0)
            .map_err(|_| QuasarError::Sumcheck("selector"))?;
        let e1 = e1.as_slice().to_vec();
        let (e1_hi, e1_lo) = e1.split_at(params.log_index());

        // The terminal claim: h̃(e1) · (ũ(e1) + σ·q̃(e1)), with the block
        // variables of ũ, q̃ resolved by the block-selector weights.
        let selector_at_e1: F = indices
            .iter()
            .enumerate()
            .map(|(k, &idx)| tau_powers[k] * eq_at_index(idx, params.log_codeword_len(), &e1))
            .sum();
        let block_selector = eq_table(e1_hi);
        let combined_at_e1: F = (0..c)
            .map(|i| {
                block_selector[i]
                    * (proof.block_evals_at_e1[i] + sigma * proof.block_evals_at_e1[c + i])
            })
            .sum();
        if selector_final != selector_at_e1 * combined_at_e1 {
            return Err(QuasarError::SelectorTerminal);
        }
        challenger.observe_algebra_slice(&proof.block_evals_at_e1);

        // --- QA-encoding check. ---
        let gamma: Vec<F> = (0..n2).map(|_| challenger.sample_algebra_element()).collect();
        let alpha: F = challenger.sample_algebra_element();
        let beta: F = challenger.sample_algebra_element();
        challenger.observe_algebra_slice(&proof.parity_evals_at_gamma);

        if proof.encoding_outer.claimed_sum
            != batched_parity_claim(c, alpha, beta, &proof.parity_evals_at_gamma)
        {
            return Err(QuasarError::EncodingClaim);
        }
        let (e2, outer_final) = proof
            .encoding_outer
            .verify(challenger, n2, 3, 0)
            .map_err(|_| QuasarError::Sumcheck("encoding outer"))?;
        let e2 = e2.as_slice().to_vec();

        // H̃ is evaluated locally in O(n₂) — the whole point of the WHT
        // formulation. The coefficient evaluations come from the preprocessed
        // commitment, and the witness evaluation is the inner sumcheck's claim.
        let coeff_at_e2: F = powers(alpha, params.num_coeff_vectors())
            .iter()
            .zip(&proof.coeff_evals_at_e2)
            .map(|(&a, &e)| a * e)
            .sum();
        if outer_final
            != wht_mle_eval(&gamma, &e2) * coeff_at_e2 * proof.encoding_inner.claimed_sum
        {
            return Err(QuasarError::EncodingOuterTerminal);
        }
        challenger.observe_algebra_slice(&proof.coeff_evals_at_e2);

        let (e3, inner_final) = proof
            .encoding_inner
            .verify(challenger, n2, 2, 0)
            .map_err(|_| QuasarError::Sumcheck("encoding inner"))?;
        let e3 = e3.as_slice().to_vec();
        if inner_final
            != wht_mle_eval(&e2, &e3)
                * (proof.input_evals_at_e3[0] + beta * proof.input_evals_at_e3[1])
        {
            return Err(QuasarError::EncodingInnerTerminal);
        }
        challenger.observe_algebra_slice(&proof.input_evals_at_e3);

        // COST ABLATION ONLY — stops after the QA-encoding sumchecks, before
        // claim reduction and the two batched Basefold openings.
        #[cfg(feature = "ablate-after-encoding")]
        return Ok(());

        // --- Claim reduction and the two batched openings. ---
        let claims = self.collect_claims(
            e1_lo,
            &gamma,
            &e2,
            &e3,
            z_lo,
            value,
            &proof.block_evals_at_e1,
            &proof.parity_evals_at_gamma,
            &proof.input_evals_at_e3,
            &proof.coeff_evals_at_e2,
        );

        let rho: F = challenger.sample_algebra_element();
        if proof.reduction.claimed_sum != batched_claim_value(&claims, rho) {
            return Err(QuasarError::ReductionClaim);
        }
        let (zeta, reduction_final) = proof
            .reduction
            .verify(challenger, n2, 2, 0)
            .map_err(|_| QuasarError::Sumcheck("reduction"))?;
        let zeta = zeta.as_slice().to_vec();

        let rho_powers = powers(rho, claims.len());
        let expected_reduction: F = claims
            .iter()
            .zip(&rho_powers)
            .map(|(claim, &w)| {
                w * eq_eval(&claim.point, &zeta) * proof.evals_at_zeta[self.flat_index(claim.target)]
            })
            .sum();
        if reduction_final != expected_reduction {
            return Err(QuasarError::ReductionTerminal);
        }
        challenger.observe_algebra_slice(&proof.evals_at_zeta);

        // COST ABLATION ONLY — stops after claim reduction, before the two
        // batched Basefold openings.
        #[cfg(feature = "ablate-after-reduction")]
        return Ok(());

        #[allow(unreachable_code)]
        let eta: F = challenger.sample_algebra_element();
        let online_etas = powers(eta, params.num_online_blocks());
        let coeff_etas = powers(eta, params.num_coeff_vectors());

        #[cfg(feature = "count-ops")]
        let probe_before = crate::opening_probe::merkle_snapshot();
        let online_value = self
            .basefold
            .verify_batch_opening(
                &proof.blocks_commitment,
                &online_etas,
                &zeta,
                &proof.online_opening,
                challenger,
            )
            .map_err(QuasarError::Basefold)?;
        #[cfg(feature = "count-ops")]
        let probe_mid = crate::opening_probe::merkle_snapshot();
        #[cfg(feature = "count-ops")]
        crate::opening_probe::record_online(probe_before, probe_mid);

        // COST ABLATION ONLY — stops after the FIRST (online) Basefold opening.
        #[cfg(feature = "ablate-after-online-opening")]
        return Ok(());

        #[allow(unreachable_code)]
        let coeff_value = self
            .basefold
            .verify_batch_opening(
                &self.coeff_commitment,
                &coeff_etas,
                &zeta,
                &proof.coeff_opening,
                challenger,
            )
            .map_err(QuasarError::Basefold)?;

        #[cfg(feature = "count-ops")]
        crate::opening_probe::record_coeff(probe_mid, crate::opening_probe::merkle_snapshot());

        // COST ABLATION ONLY — stops after BOTH Basefold openings, before the
        // final batched-value equality checks.
        #[cfg(feature = "ablate-after-both-openings")]
        {
            let _ = (online_value, coeff_value);
            return Ok(());
        }

        let expected_online: F = online_etas
            .iter()
            .zip(&proof.evals_at_zeta[..params.num_online_blocks()])
            .map(|(&e, &y)| e * y)
            .sum();
        let expected_coeff: F = coeff_etas
            .iter()
            .zip(&proof.evals_at_zeta[params.num_online_blocks()..])
            .map(|(&e, &y)| e * y)
            .sum();
        if online_value != expected_online || coeff_value != expected_coeff {
            return Err(QuasarError::OpeningValue);
        }

        Ok(())
    }

    // ---- shared helpers (identical on both sides) ----

    /// `out[j] = Σ_i weights[i]·M[i][j]`, from the transposed message matrix.
    ///
    /// `O(N)` and the dominant non-hashing cost of the evaluation phase; each
    /// output column is an independent contiguous dot product.
    fn combine_rows(&self, msg: &RowMajorMatrix<F>, weights: &[F]) -> Vec<F> {
        let rows = self.params.num_rows();
        msg.values
            .par_chunks_exact(rows)
            .map(|row| dot(row, weights))
            .collect()
    }

    /// `Σ_i α^i · (2^{-n₂}a_{i+2}')`, the batched public coefficient vector.
    fn combine_coeffs(&self, alpha: F) -> Vec<F> {
        let cols = self.params.num_cols();
        let mut out = vec![F::ZERO; cols];
        for (coeff, w) in self
            .code()
            .scaled_coeffs()
            .iter()
            .zip(powers(alpha, self.params.num_coeff_vectors()))
        {
            out.par_iter_mut()
                .zip(coeff.par_iter())
                .for_each(|(o, &e)| *o += w * e);
        }
        out
    }

    /// The sparse selector table: `τ^k` at each sampled column index.
    ///
    /// Repeated indices accumulate, which is exactly what makes the right-hand
    /// side `Σ_k τ^k⟨w, C(·,i_k)⟩` correct for a multiset of samples.
    fn selector_table(&self, indices: &[usize], tau: F) -> Vec<F> {
        let mut table = vec![F::ZERO; self.params.codeword_len()];
        for (weight, &idx) in powers(tau, indices.len()).iter().zip(indices) {
            table[idx] += *weight;
        }
        table
    }

    /// Sample `num_queries` column indices from the transcript.
    fn sample_columns<Ch: CanSampleBits<usize>>(&self, challenger: &mut Ch) -> Vec<usize> {
        let n = self.params.codeword_len();
        let bits = log2_ceil_usize(n);
        let mut out = Vec::with_capacity(self.params.num_queries);
        while out.len() < self.params.num_queries {
            let idx = challenger.sample_bits(bits);
            if idx < n {
                out.push(idx);
            }
        }
        out
    }

    /// Position of a target's evaluation in `evals_at_zeta`.
    const fn flat_index(&self, target: Target) -> usize {
        match target {
            Target::Online(j) => j,
            Target::Coeff(i) => self.params.num_online_blocks() + i,
        }
    }

    /// Every evaluation claim the batched openings must discharge, in a fixed
    /// order both sides reproduce.
    #[allow(clippy::too_many_arguments)]
    fn collect_claims(
        &self,
        e1_lo: &[F],
        gamma: &[F],
        e2: &[F],
        e3: &[F],
        z_lo: &[F],
        value: F,
        block_evals_at_e1: &[F],
        parity_evals_at_gamma: &[F],
        input_evals_at_e3: &[F; 2],
        coeff_evals_at_e2: &[F],
    ) -> Vec<Claim<F>> {
        let c = self.params.index;
        let mut claims = Vec::new();

        // Selector terminal: every online block at e1_lo.
        for (j, &value) in block_evals_at_e1.iter().enumerate() {
            claims.push(Claim {
                target: Target::Online(j),
                point: e1_lo.to_vec(),
                value,
            });
        }
        // Encoding claim: the parity blocks at gamma.
        for i in 1..c {
            claims.push(Claim {
                target: Target::Online(i),
                point: gamma.to_vec(),
                value: parity_evals_at_gamma[i - 1],
            });
            claims.push(Claim {
                target: Target::Online(c + i),
                point: gamma.to_vec(),
                value: parity_evals_at_gamma[c - 2 + i],
            });
        }
        // Inner terminal: the two input blocks at e3.
        claims.push(Claim {
            target: Target::Online(0),
            point: e3.to_vec(),
            value: input_evals_at_e3[0],
        });
        claims.push(Claim {
            target: Target::Online(c),
            point: e3.to_vec(),
            value: input_evals_at_e3[1],
        });
        // The evaluation being proved: q₁(z_lo) = y, free because the code is
        // systematic (q₁ = p = z_leftᵗM).
        claims.push(Claim {
            target: Target::Online(c),
            point: z_lo.to_vec(),
            value,
        });
        // Preprocessed coefficients at the outer terminal point.
        for (i, &value) in coeff_evals_at_e2.iter().enumerate() {
            claims.push(Claim {
                target: Target::Coeff(i),
                point: e2.to_vec(),
                value,
            });
        }
        claims
    }

    /// Per-vector weight tables `W_j = Σ_{claims on j} ρ^k·eq(point_k, ·)`.
    fn claim_weights(&self, claims: &[Claim<F>], rho: F) -> ClaimWeights<F> {
        let cols = self.params.num_cols();
        let mut weights = ClaimWeights {
            online: vec![vec![F::ZERO; cols]; self.params.num_online_blocks()],
            coeff: vec![vec![F::ZERO; cols]; self.params.num_coeff_vectors()],
        };
        for (claim, w) in claims.iter().zip(powers(rho, claims.len())) {
            let table = eq_table(&claim.point);
            let dst = match claim.target {
                Target::Online(j) => &mut weights.online[j],
                Target::Coeff(i) => &mut weights.coeff[i],
            };
            dst.par_iter_mut()
                .zip(table.par_iter())
                .for_each(|(d, &e)| *d += w * e);
        }
        weights
    }
}

/// Grouped claim weights, one table per committed vector.
struct ClaimWeights<F> {
    online: Vec<Vec<F>>,
    coeff: Vec<Vec<F>>,
}

/// `Σ_i α^i·(u_{i+1}(γ) + β·q_{i+1}(γ))` — the outer sumcheck's claimed sum.
fn batched_parity_claim<F: Field>(index: usize, alpha: F, beta: F, parity_evals: &[F]) -> F {
    let c = index;
    powers(alpha, c - 1)
        .iter()
        .enumerate()
        .map(|(i, &w)| w * (parity_evals[i] + beta * parity_evals[c - 1 + i]))
        .sum()
}

/// `Σ_k ρ^k·y_k` — the claim-reduction sumcheck's claimed sum.
fn batched_claim_value<F: Field>(claims: &[Claim<F>], rho: F) -> F {
    claims
        .iter()
        .zip(powers(rho, claims.len()))
        .map(|(claim, w)| w * claim.value)
        .sum()
}
