//! A standalone tensor (Ligero/Brakedown) polynomial commitment for a single
//! multilinear polynomial (M3a).
//!
//! Field-agnostic: `Val: Field`, challenges in `EF: ExtensionField<Val>`, a
//! linear code over `Val`, and any Plonky3 [`Mmcs`] over `Val`. No FFT / no
//! two-adicity.
//!
//! ## Construction
//!
//! A polynomial `p` over `m` variables has `2^m` evaluations. Reshape them into
//! a `rows × cols` matrix `M` (`rows = 2^log_rows`, `cols = code.message_len()`,
//! `rows·cols = 2^m`), where `M[i][j] = evals[i·cols + j]`.
//!
//! * **Commit**: encode each row of `M` with the code → `M̂` (`rows × n`); Merkle-
//!   commit the *columns* of `M̂` (we commit `M̂ᵀ`, so leaf `t` = column `t`).
//! * **Open** at point `z = (z_hi ‖ z_lo)` (`log_rows` high vars, the rest low):
//!   with `b = eq(·, z_hi)` (len `rows`) and `a = eq(·, z_lo)` (len `cols`),
//!   `p(z) = ⟨b, M, a⟩`. Send the combined row `w = bᵀM` (len `cols`), a random
//!   proximity row `w' = γᵀM`, and `t` queried columns with Merkle proofs.
//! * **Verify**: for each queried column `c`, check `⟨b, M̂[:,c]⟩ = Encode(w)[c]`
//!   (evaluation) and `⟨γ, M̂[:,c]⟩ = Encode(w')[c]` (proximity / well-formedness),
//!   and that `⟨w, a⟩` equals the claimed value.
//!
//! ## Scope (M3a)
//!
//! This proves a *plain point* evaluation `p(z) = v`. Integrating with
//! `multi-stark` (the non-tensor "next" view) needs a sumcheck-batching layer on
//! top (M3b). Soundness parameters (query count vs. code distance) are M4;
//! `num_queries` here is a caller-supplied knob.
//!
//! ## Soundness note: `num_queries` vs. code distance
//!
//! For a linear code with minimum *relative* distance `δ`, the column-query test
//! is a tensor/Ligero-style proximity + consistency argument. Roughly, a cheating
//! prover whose committed matrix is `δ'`-far from a valid codeword survives a
//! single random column query with probability at most `1 − δ'` (the query lands
//! in an agreeing position), so `num_queries` independent columns give a
//! consistency soundness error on the order of `(1 − δ)^num_queries` (with the
//! test's proximity gap capping the usable `δ`). On top of that sits an additive
//! **proximity-test term** — the probability that the random linear combination
//! `w' = γᵀM` masks a far-from-codeword row — on the order of the field-size
//! reciprocal / the proximity gap (see the Ligero and Brakedown analyses).
//!
//! Concretely: `soundness_error ≈ (1 − δ)^num_queries + proximity_term`. Choosing
//! `num_queries` for a target security level (e.g. 100 bits) therefore requires
//! the code's *proven* relative distance `δ`. For the Reed–Solomon base
//! (`δ = 1 − k/n`) this is exact; for the recursive Brakedown code the proven `δ`
//! needs the paper's exact degree formulas (see
//! [`BrakedownParams`](crate::BrakedownParams)). Deriving `num_queries` from a
//! target bit-security is therefore **deferred / left configurable**;
//! `num_queries` is a caller-supplied knob here.

use alloc::vec::Vec;
use core::marker::PhantomData;

use p3_challenger::{CanObserve, CanSampleBits, FieldChallenger};
use p3_commit::Mmcs;
use p3_field::{ExtensionField, Field};
use p3_matrix::Dimensions;
use p3_matrix::dense::RowMajorMatrix;
use p3_multilinear_util::point::Point;
use p3_util::log2_ceil_usize;
use serde::{Deserialize, Serialize};

use crate::LinearCodeExt;

/// A tensor-code multilinear PCS over base field `Val`, challenge field `EF`,
/// linear code `C`, and Merkle scheme `M`.
#[derive(Clone)]
pub struct TensorPcs<Val, EF, C, M> {
    code: C,
    mmcs: M,
    /// Number of high (row) variables; `rows = 2^log_rows`.
    log_rows: usize,
    /// Number of columns queried in the proximity + evaluation test.
    num_queries: usize,
    _phantom: PhantomData<(Val, EF)>,
}

/// Prover state retained between commit and open.
pub struct TensorProverData<Val: Send + Sync + Clone, M: Mmcs<Val>> {
    /// The message matrix in `cols × rows` layout (`msg[j][i] = M[i][j]`), used
    /// to form the combined rows.
    msg: RowMajorMatrix<Val>,
    /// Merkle prover data over `M̂ᵀ` (the committed encoded columns).
    mmcs_data: M::ProverData<RowMajorMatrix<Val>>,
}

// Hand-written so the `Clone` bound lands on the Merkle prover data associated
// type rather than on `M` itself (a derive cannot express that).
impl<Val: Send + Sync + Clone, M: Mmcs<Val>> Clone for TensorProverData<Val, M>
where
    M::ProverData<RowMajorMatrix<Val>>: Clone,
{
    fn clone(&self) -> Self {
        Self {
            msg: self.msg.clone(),
            mmcs_data: self.mmcs_data.clone(),
        }
    }
}

/// An opening proof for a single point evaluation.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(
    serialize = "Val: Serialize, EF: Serialize, M::Proof: Serialize",
    deserialize = "Val: serde::de::DeserializeOwned, EF: serde::de::DeserializeOwned, M::Proof: serde::de::DeserializeOwned"
))]
pub struct TensorProof<Val: Send + Sync + Clone, EF, M: Mmcs<Val>> {
    /// The single combined row `w = bᵀM` (length `cols`). With a **random** opening point, `b`
    /// (the eq-weights of the point's row variables) is a Diamond–Posen tensor combination that
    /// serves as BOTH the evaluation and the proximity check — no separate proximity row is needed.
    pub w: Vec<EF>,
    /// The `t` queried columns of `M̂`: `openings[k][0]` is column `indices[k]`,
    /// all under ONE pruned authentication path. The queries share most of their
    /// ancestors, so the shared siblings are recomputable and need not be sent.
    pub openings: Vec<Vec<Vec<Val>>>,
    /// The pruned paths covering every queried column.
    pub opening_proof: M::MultiProof,
}

/// Verification failures.
#[derive(Debug)]
pub enum TensorPcsError<E> {
    /// A Merkle column opening failed.
    Merkle(E),
    /// The combined consistency+proximity check `⟨b, col⟩ = Encode(w)[c]` failed at a query.
    EvalCheckFailed { query: usize },
    /// The claimed value did not equal `⟨w, a⟩`.
    ClaimedValueMismatch,
    /// The proof carried the wrong number of column openings.
    WrongOpeningCount { expected: usize, got: usize },
}

impl<Val, EF, C, M> TensorPcs<Val, EF, C, M>
where
    Val: Field,
    EF: ExtensionField<Val>,
    C: LinearCodeExt<Val, EF>,
    M: Mmcs<Val>,
{
    /// Build a tensor PCS. `log_rows` sets the row/column split; the polynomial
    /// committed must have `2^log_rows · code.message_len()` evaluations.
    pub fn new(code: C, mmcs: M, log_rows: usize, num_queries: usize) -> Self {
        Self {
            code,
            mmcs,
            log_rows,
            num_queries,
            _phantom: PhantomData,
        }
    }

    fn rows(&self) -> usize {
        1 << self.log_rows
    }

    fn cols(&self) -> usize {
        self.code.message_len()
    }

    fn codeword_len(&self) -> usize {
        self.code.codeword_len()
    }

    /// Number of variables of the committed multilinear polynomial.
    pub fn num_variables(&self) -> usize {
        self.log_rows + log2_ceil_usize(self.cols())
    }

    /// Commit to a multilinear polynomial given by its `2^m` evaluations.
    pub fn commit<Ch>(
        &self,
        evals: Vec<Val>,
        challenger: &mut Ch,
    ) -> (M::Commitment, TensorProverData<Val, M>)
    where
        Ch: FieldChallenger<Val> + CanObserve<M::Commitment>,
    {
        let (rows, cols) = (self.rows(), self.cols());
        assert_eq!(evals.len(), rows * cols, "eval count must be rows*cols");

        // Natural reshape M (rows × cols), then transpose to msg (cols × rows)
        // so the code (which encodes along its height) encodes each M-row.
        let m_natural = RowMajorMatrix::new(evals, cols);
        let msg = m_natural.transpose(); // cols × rows
        // n × rows  (== M̂ᵀ; leaf t = column t of M̂). Spans surface the two
        // commit stages — linear-code encoding vs Merkle build — in a prover's
        // tracing tree (e.g. `eddsa-example`), where they nest under `prove`.
        let mhat_t = tracing::info_span!("encode").in_scope(|| self.code.encode(&msg));

        let (commitment, mmcs_data) =
            tracing::info_span!("merkle_commit").in_scope(|| self.mmcs.commit_matrix(mhat_t));
        challenger.observe(commitment.clone());
        (commitment, TensorProverData { msg, mmcs_data })
    }

    /// Open the committed polynomial at `point`, returning the evaluation and a proof.
    pub fn open<Ch>(
        &self,
        prover_data: &TensorProverData<Val, M>,
        point: &Point<EF>,
        challenger: &mut Ch,
    ) -> (EF, TensorProof<Val, EF, M>)
    where
        Ch: FieldChallenger<Val> + CanSampleBits<usize>,
    {
        assert_eq!(point.num_variables(), self.num_variables(), "point arity");

        let (b, a) = self.eq_weights(point);

        // Consolidated (Diamond–Posen): a single combined row w = bᵀM. Because the opening point
        // is random, b (its row-variable eq-weights, a tensor combination) doubles as the proximity
        // combination — the one per-column check below is both consistency AND proximity.
        let w = self.combine_rows(&prover_data.msg, &b);
        challenger.observe_algebra_slice(&w);

        // Claimed evaluation p(z) = ⟨w, a⟩.
        let value = dot(&w, &a);

        // Sample query columns and open them.
        let indices = self.sample_queries(challenger);
        let (openings, opening_proof) = self
            .mmcs
            .open_multi_batch(&indices, &prover_data.mmcs_data);

        (value, TensorProof { w, openings, opening_proof })
    }

    /// Verify an opening proof for `point`, checking the claimed `value`.
    pub fn verify<Ch>(
        &self,
        commitment: &M::Commitment,
        point: &Point<EF>,
        value: EF,
        proof: &TensorProof<Val, EF, M>,
        challenger: &mut Ch,
    ) -> Result<(), TensorPcsError<M::Error>>
    where
        Ch: FieldChallenger<Val> + CanObserve<M::Commitment> + CanSampleBits<usize>,
    {
        assert_eq!(point.num_variables(), self.num_variables(), "point arity");
        let (rows, cols) = (self.rows(), self.cols());
        let n = self.codeword_len();

        // Replay the transcript exactly as the prover drove it.
        challenger.observe(commitment.clone());
        challenger.observe_algebra_slice(&proof.w);
        let indices = self.sample_queries(challenger);

        if proof.openings.len() != indices.len() {
            return Err(TensorPcsError::WrongOpeningCount {
                expected: indices.len(),
                got: proof.openings.len(),
            });
        }
        assert_eq!(proof.w.len(), cols);

        let (b, a) = self.eq_weights(point);

        // Claimed value must equal ⟨w, a⟩ (ties `value` to the combined row `w`).
        if value != dot(&proof.w, &a) {
            return Err(TensorPcsError::ClaimedValueMismatch);
        }

        // Encode the combined row once; symbol `c` is Encode(w)[c].
        let enc_w = self.encode_row(&proof.w);

        // The committed matrix M̂ᵀ has n rows (leaves), each of width `rows`.
        let dims = [Dimensions {
            width: rows,
            height: n,
        }];

        if proof.openings.iter().any(|v| v.len() != 1 || v[0].len() != rows) {
            return Err(TensorPcsError::WrongOpeningCount {
                expected: indices.len(),
                got: proof.openings.len(),
            });
        }
        // One pruned multi-opening covers every queried column.
        self.mmcs
            .verify_multi_batch(commitment, &dims, &indices, &proof.openings, &proof.opening_proof)
            .map_err(TensorPcsError::Merkle)?;

        for (query, (&idx, opening)) in indices.iter().zip(&proof.openings).enumerate() {
            // leaf = column `idx` of M̂: consolidated consistency + proximity,
            // ⟨b, leaf⟩ == Encode(w)[idx].
            if dot_base(&b, &opening[0]) != enc_w[idx] {
                return Err(TensorPcsError::EvalCheckFailed { query });
            }
        }
        Ok(())
    }

    /// Verify an opening proof for `point` **without** re-absorbing the
    /// commitment, returning the value `v = ⟨w, a⟩` the proof binds `point` to.
    ///
    /// This is the integration-facing variant of [`verify`](Self::verify): the
    /// caller (an outer multilinear-PCS adapter) has already absorbed the
    /// commitment into the shared transcript, and checks the returned value
    /// against its own batched claim instead of passing a pre-committed `value`.
    ///
    /// It runs the same eta / w / w' / query-index transcript sequence as
    /// [`open`](Self::open) and the same Merkle + evaluation + proximity checks
    /// as [`verify`](Self::verify), but skips the commitment observe and the
    /// `⟨w, a⟩ == value` equality (returning `⟨w, a⟩` for the caller to check).
    pub fn verify_opening<Ch>(
        &self,
        commitment: &M::Commitment,
        point: &Point<EF>,
        proof: &TensorProof<Val, EF, M>,
        challenger: &mut Ch,
    ) -> Result<EF, TensorPcsError<M::Error>>
    where
        Ch: FieldChallenger<Val> + CanSampleBits<usize>,
    {
        assert_eq!(point.num_variables(), self.num_variables(), "point arity");
        let (rows, cols) = (self.rows(), self.cols());
        let n = self.codeword_len();

        // Replay the post-commitment transcript exactly as the prover drove it in
        // `open` (the commitment itself was absorbed by the outer protocol).
        challenger.observe_algebra_slice(&proof.w);
        let indices = self.sample_queries(challenger);

        if proof.openings.len() != indices.len() {
            return Err(TensorPcsError::WrongOpeningCount {
                expected: indices.len(),
                got: proof.openings.len(),
            });
        }
        assert_eq!(proof.w.len(), cols);

        let (b, a) = self.eq_weights(point);

        // Encode the combined row once; symbol `c` is Encode(w)[c].
        let enc_w = self.encode_row(&proof.w);

        let dims = [Dimensions {
            width: rows,
            height: n,
        }];

        if proof.openings.iter().any(|v| v.len() != 1 || v[0].len() != rows) {
            return Err(TensorPcsError::WrongOpeningCount {
                expected: indices.len(),
                got: proof.openings.len(),
            });
        }
        // One pruned multi-opening covers every queried column.
        self.mmcs
            .verify_multi_batch(commitment, &dims, &indices, &proof.openings, &proof.opening_proof)
            .map_err(TensorPcsError::Merkle)?;

        for (query, (&idx, opening)) in indices.iter().zip(&proof.openings).enumerate() {
            // leaf = column `idx` of M̂: consolidated consistency + proximity,
            // ⟨b, leaf⟩ == Encode(w)[idx].
            if dot_base(&b, &opening[0]) != enc_w[idx] {
                return Err(TensorPcsError::EvalCheckFailed { query });
            }
        }

        // The proof binds `point` to this value; the caller cross-checks it.
        Ok(dot(&proof.w, &a))
    }

    // ---- helpers ----

    /// Split `point` into row/column equality-weight tables `(b, a)`.
    fn eq_weights(&self, point: &Point<EF>) -> (Vec<EF>, Vec<EF>) {
        let coords = point.as_slice();
        let (z_hi, z_lo) = coords.split_at(self.log_rows);
        (eq_table(z_hi), eq_table(z_lo))
    }

    /// Combined row: `out[j] = Σ_i weights[i] · msg[j][i]` (msg is cols × rows).
    fn combine_rows(&self, msg: &RowMajorMatrix<Val>, weights: &[EF]) -> Vec<EF> {
        let (rows, cols) = (self.rows(), self.cols());
        (0..cols)
            .map(|j| {
                let row = &msg.values[j * rows..(j + 1) * rows];
                let mut acc = EF::ZERO;
                for (i, &m) in row.iter().enumerate() {
                    acc += weights[i] * m;
                }
                acc
            })
            .collect()
    }

    /// Encode an `EF` combined row (length cols) to its length-`n` codeword.
    fn encode_row(&self, row: &[EF]) -> Vec<EF> {
        let m = RowMajorMatrix::new(row.to_vec(), 1);
        self.code.encode_ext(&m).values
    }

    /// Sample `num_queries` column indices in `[0, n)` from the transcript.
    fn sample_queries<Ch: CanSampleBits<usize>>(&self, challenger: &mut Ch) -> Vec<usize> {
        let n = self.codeword_len();
        let bits = log2_ceil_usize(n);
        let mut out = Vec::with_capacity(self.num_queries);
        while out.len() < self.num_queries {
            let idx = challenger.sample_bits(bits);
            if idx < n {
                out.push(idx);
            }
        }
        out
    }
}

/// The equality-weight table `eq(x, z)` for all `x ∈ {0,1}^{|z|}`, big-endian
/// (coordinate 0 is the most significant bit of `x`).
fn eq_table<EF: Field>(z: &[EF]) -> Vec<EF> {
    let mut table = alloc::vec![EF::ONE];
    for &zi in z {
        let mut next = Vec::with_capacity(table.len() * 2);
        for &t in &table {
            next.push(t * (EF::ONE - zi)); // bit = 0
            next.push(t * zi); // bit = 1
        }
        table = next;
    }
    table
}

/// `Σ_i x[i] · y[i]` over `EF`.
fn dot<EF: Field>(x: &[EF], y: &[EF]) -> EF {
    x.iter().zip(y).map(|(&a, &b)| a * b).sum()
}

/// `Σ_i x[i] · y[i]` where `x ∈ EF`, `y ∈ Val`.
fn dot_base<Val: Field, EF: ExtensionField<Val>>(x: &[EF], y: &[Val]) -> EF {
    x.iter().zip(y).map(|(&a, &b)| a * b).sum()
}
