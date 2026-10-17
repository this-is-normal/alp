//! M2-auth: the **succinct** Basefold PCS — Merkle-committed folded oracles + a `q_bf`-query phase
//! driven by a real Fiat–Shamir challenger. This is the polylog-verifier version (the M2-core in
//! `pcs.rs` sent full oracles; here the verifier only checks `q_bf` query paths).
//!
//! Single field (`Challenge = Val`), matching alp's eddsa config (`EF = Ed25519 = F`).

// `vec!` is only used by the tests below; the library builds without it.
#[cfg(test)]
use alloc::vec;
use alloc::vec::Vec;

use p3_challenger::{CanObserve, CanSampleBits, FieldChallenger, GrindingChallenger};
use p3_commit::Mmcs;
use p3_field::Field;
use p3_matrix::Dimensions;
use p3_matrix::dense::RowMajorMatrix;
use p3_maybe_rayon::prelude::*;
use p3_util::log2_ceil_usize;
use serde::{Deserialize, Serialize};

use crate::foldable_code::{FoldableCode, RandomFoldableCode};
use crate::pcs::{bind_msb, eq_table, eval_eq, interp_deg2};

/// Prover data retained between `commit` and `open`.
pub struct BasefoldProverData<F: Field, M: Mmcs<F>> {
    evals: Vec<F>,
    o0_cw: Vec<F>,
    o0_data: M::ProverData<RowMajorMatrix<F>>,
}

impl<F: Field, M: Mmcs<F>> BasefoldProverData<F, M> {
    /// The committed stacked-polynomial evaluations (needed by the multi-stark adapter's sumcheck).
    pub fn evals_slice(&self) -> &[F] {
        &self.evals
    }
}

impl<F: Field, M: Mmcs<F>> Clone for BasefoldProverData<F, M>
where
    M::ProverData<RowMajorMatrix<F>>: Clone,
{
    fn clone(&self) -> Self {
        Self { evals: self.evals.clone(), o0_cw: self.o0_cw.clone(), o0_data: self.o0_data.clone() }
    }
}

/// A succinct Basefold evaluation proof.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(
    serialize = "F: Serialize, M::Commitment: Serialize, M::Proof: Serialize",
    deserialize = "F: serde::de::DeserializeOwned, M::Commitment: serde::de::DeserializeOwned, M::Proof: serde::de::DeserializeOwned"
))]
pub struct BasefoldProof<F: Field, M: Mmcs<F>> {
    /// Sumcheck univariates `h_0..h_{d-1}` at points 0,1,2.
    pub sumcheck: Vec<[F; 3]>,
    /// Merkle roots of the folded oracles `O_1..O_{d-1}` (`O_0` is the PCS commitment).
    pub fold_roots: Vec<M::Commitment>,
    /// The final folded value `v = f(r)`.
    pub final_value: F,
    /// Proof-of-work witness bound just before the query indices were drawn.
    /// `F::ZERO` when `pow_bits == 0`.
    pub pow_witness: F,
    /// Per round `i` (0..d): every sampled folding-pair of `O_i`, under ONE pruned
    /// authentication path.
    ///
    /// Pruning matters here more than anywhere else in the scheme: the `q` queries
    /// into a round-`i` oracle share most of their ancestors, and the shared
    /// siblings are recomputable, so sending them is pure waste. Unpruned, these
    /// paths are the bulk of a Quasar proof.
    pub queries: Vec<RoundOpening<F, M>>,
}

/// Every query's opened row for one folded oracle, with a single shared proof.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(
    serialize = "F: Serialize, M::MultiProof: Serialize",
    deserialize = "F: serde::de::DeserializeOwned, M::MultiProof: serde::de::DeserializeOwned"
))]
pub struct RoundOpening<F: Field, M: Mmcs<F>> {
    /// `opened_values[q][0]` is the folding pair (or batched pairs) at query `q`.
    pub opened_values: Vec<Vec<Vec<F>>>,
    /// The pruned authentication paths covering all of them.
    pub proof: M::MultiProof,
}

/// Additive accumulator for a multi-value reduction.
///
/// `p3-maybe-rayon`'s two backends disagree on `reduce`: rayon's takes an
/// identity plus a combiner, `core::iter`'s takes only a combiner. `sum` has the
/// same shape in both, so accumulating through this type is what lets one
/// expression compile serial and parallel.
#[derive(Copy, Clone)]
struct Sums<F, const N: usize>([F; N]);

impl<F: Field, const N: usize> Default for Sums<F, N> {
    fn default() -> Self {
        Self([F::ZERO; N])
    }
}

impl<F: Field, const N: usize> core::ops::Add for Sums<F, N> {
    type Output = Self;

    fn add(mut self, other: Self) -> Self {
        for (a, b) in self.0.iter_mut().zip(other.0) {
            *a += b;
        }
        self
    }
}

impl<F: Field, const N: usize> core::iter::Sum for Sums<F, N> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), |a, b| a + b)
    }
}

/// Reshape a codeword `cw` (length `n`) into an `(n/2) × 2` matrix whose row `j` is the folding
/// pair `[cw[j], cw[j+n/2]]`, so one Merkle path opens both siblings of a fold.
fn pair_matrix<F: Field>(cw: &[F]) -> RowMajorMatrix<F> {
    pair_matrix_many(core::slice::from_ref(&cw))
}

/// The batched form: `k` equal-length codewords side by side, as an `(n/2) × 2k`
/// matrix whose row `j` is `[cw_0[j], cw_0[j+n/2], cw_1[j], cw_1[j+n/2], …]`.
///
/// One Merkle path therefore opens the folding pair of *every* batched codeword
/// at that position, which is what makes a batched opening cost one path per
/// query rather than `k`.
fn pair_matrix_many<F: Field>(cws: &[&[F]]) -> RowMajorMatrix<F> {
    assert!(!cws.is_empty(), "need at least one codeword");
    let n = cws[0].len();
    assert!(cws.iter().all(|cw| cw.len() == n), "codeword lengths differ");
    let half = n / 2;
    let width = 2 * cws.len();
    let mut flat = F::zero_vec(n * cws.len());
    flat.par_chunks_exact_mut(width)
        .enumerate()
        .for_each(|(j, row)| {
            for (k, cw) in cws.iter().enumerate() {
                row[2 * k] = cw[j];
                row[2 * k + 1] = cw[j + half];
            }
        });
    RowMajorMatrix::new(flat, width)
}

/// What the shared folding engine hands back: the final folded value, the
/// per-round sumcheck univariates, and the folded oracles' roots and prover
/// data.
type FoldingOutput<F, M> = (
    F,
    Vec<[F; 3]>,
    Vec<<M as Mmcs<F>>::Commitment>,
    Vec<<M as Mmcs<F>>::ProverData<RowMajorMatrix<F>>>,
);

/// Prover data for a batch commitment to `k` equal-length polynomials under one
/// Merkle root — the shape a batched opening needs (see
/// [`BasefoldPcs::commit_many`]).
pub struct BasefoldBatchProverData<F: Field, M: Mmcs<F>> {
    polys: Vec<Vec<F>>,
    codewords: Vec<Vec<F>>,
    data: M::ProverData<RowMajorMatrix<F>>,
}

impl<F: Field, M: Mmcs<F>> BasefoldBatchProverData<F, M> {
    /// Number of batched polynomials.
    pub const fn len(&self) -> usize {
        self.polys.len()
    }

    /// Whether the batch is empty (never, in practice: `commit_many` rejects it).
    pub const fn is_empty(&self) -> bool {
        self.polys.is_empty()
    }

    /// The `j`-th committed evaluation vector.
    pub fn poly(&self, j: usize) -> &[F] {
        &self.polys[j]
    }
}

/// The succinct Basefold PCS over a foldable code `code` and a Merkle commitment scheme `mmcs`.
pub struct BasefoldPcs<F: Field, M: Mmcs<F>> {
    pub code: RandomFoldableCode<F>,
    pub mmcs: M,
    pub num_queries: usize,
    /// Proof-of-work bits ground immediately before the query indices are
    /// sampled. Each bit buys one bit of query soundness, so `num_queries` can
    /// drop by `pow_bits / −log₂(1−e)` for the same target.
    ///
    /// Grinding *after* the final value is absorbed and *before* the indices are
    /// drawn is what makes it bind: a prover who dislikes its query set must
    /// redo the work to get another.
    pub pow_bits: usize,
}

impl<F: Field, M: Mmcs<F>> BasefoldPcs<F, M> {
    /// Without grinding.
    pub const fn new(code: RandomFoldableCode<F>, mmcs: M, num_queries: usize) -> Self {
        Self { code, mmcs, num_queries, pow_bits: 0 }
    }

    /// With `pow_bits` of proof of work before the query phase.
    pub const fn new_with_grinding(
        code: RandomFoldableCode<F>,
        mmcs: M,
        num_queries: usize,
        pow_bits: usize,
    ) -> Self {
        Self { code, mmcs, num_queries, pow_bits }
    }

    pub fn num_vars(&self) -> usize {
        self.code.log_len()
    }

    /// Commit to `evals` (length `2^d`): encode → Merkle-commit the folding-pair matrix of `O_0`.
    pub fn commit<Ch>(&self, evals: Vec<F>, challenger: &mut Ch) -> (M::Commitment, BasefoldProverData<F, M>)
    where
        Ch: CanObserve<M::Commitment>,
    {
        assert_eq!(evals.len(), 1 << self.num_vars());
        let o0_cw = self.code.encode(&evals);
        let (commitment, o0_data) = self.mmcs.commit_matrix(pair_matrix(&o0_cw));
        challenger.observe(commitment.clone());
        (commitment, BasefoldProverData { evals, o0_cw, o0_data })
    }

    /// Grind the proof-of-work that gates the query phase.
    fn grind<Ch>(&self, challenger: &mut Ch) -> F
    where
        Ch: GrindingChallenger<Witness = F> + Clone + Send + Sync,
    {
        crate::grind_parallel(challenger, self.pow_bits)
    }

    /// Check the witness the prover ground. Must run at exactly the point the
    /// prover ground it, or the query indices desync.
    fn check_grind<Ch>(&self, witness: F, challenger: &mut Ch) -> Result<(), &'static str>
    where
        Ch: GrindingChallenger<Witness = F>,
    {
        if self.pow_bits == 0 || challenger.check_witness(self.pow_bits, witness) {
            Ok(())
        } else {
            Err("invalid proof-of-work witness")
        }
    }

    fn query_indices<Ch: CanSampleBits<usize>>(&self, challenger: &mut Ch) -> Vec<usize> {
        let rows0 = self.code.codeword_len() / 2; // rows of O_0's pair-matrix
        let bits = log2_ceil_usize(rows0);
        let mut out = Vec::with_capacity(self.num_queries);
        while out.len() < self.num_queries {
            let idx = challenger.sample_bits(bits);
            if idx < rows0 {
                out.push(idx);
            }
        }
        out
    }

    /// Commit to `k` equal-length polynomials under a single Merkle root, so a
    /// later [`Self::open_batch`] proves any linear combination of them with one
    /// path per query.
    ///
    /// The Quasar code-switching layer commits all its length-`N₂` blocks this
    /// way: without it, the several evaluation claims per block would each need
    /// their own full Basefold proof.
    pub fn commit_many<Ch>(
        &self,
        polys: Vec<Vec<F>>,
        challenger: &mut Ch,
    ) -> (M::Commitment, BasefoldBatchProverData<F, M>)
    where
        Ch: CanObserve<M::Commitment>,
    {
        let (commitment, data) = self.commit_many_silent(polys);
        challenger.observe(commitment.clone());
        (commitment, data)
    }

    /// [`Self::commit_many`] without touching a transcript — for preprocessing,
    /// where the commitment becomes part of the verifier key rather than a
    /// prover message.
    pub fn commit_many_silent(
        &self,
        polys: Vec<Vec<F>>,
    ) -> (M::Commitment, BasefoldBatchProverData<F, M>) {
        assert!(!polys.is_empty(), "batch commitment needs at least one poly");
        assert!(
            polys.iter().all(|p| p.len() == 1 << self.num_vars()),
            "every batched poly must have 2^num_vars evaluations"
        );
        let codewords: Vec<Vec<F>> = tracing::info_span!("basefold_encode")
            .in_scope(|| polys.iter().map(|p| self.code.encode(p)).collect());
        let refs: Vec<&[F]> = codewords.iter().map(Vec::as_slice).collect();
        let (commitment, data) = tracing::info_span!("basefold_merkle")
            .in_scope(|| self.mmcs.commit_matrix(pair_matrix_many(&refs)));
        (commitment, BasefoldBatchProverData { polys, codewords, data })
    }

    /// Prove `f(z) = y`.
    pub fn open<Ch>(&self, pd: &BasefoldProverData<F, M>, z: &[F], challenger: &mut Ch) -> (F, BasefoldProof<F, M>)
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>
            + Clone
            + Send
            + Sync,
    {
        let (v, sumcheck, fold_roots, fold_data) =
            self.prove_folding(pd.evals.clone(), pd.o0_cw.clone(), z, challenger);
        let pow_witness = self.grind(challenger);
        let indices = self.query_indices(challenger);
        let queries = self.collect_queries(&indices, &pd.o0_data, &fold_data);
        (v, BasefoldProof { sumcheck, fold_roots, final_value: v, pow_witness, queries })
    }

    /// Prove `g(z) = y` for the virtual polynomial `g = Σ_j etas[j]·f_j` over a
    /// batch commitment, in one proof.
    ///
    /// Sound because the code is linear: the codeword of the combination is the
    /// combination of the codewords, and the verifier reconstructs the round-0
    /// folding pair of `g` from the batched pair it opens.
    ///
    /// # Returns
    ///
    /// `(g(z), proof)` — the claimed value the proof binds, read off the first
    /// sumcheck round as `h₀(0) + h₀(1)`. (Note this differs from [`Self::open`],
    /// whose first return is the *final folded* value `f(r)`.)
    pub fn open_batch<Ch>(
        &self,
        pd: &BasefoldBatchProverData<F, M>,
        etas: &[F],
        z: &[F],
        challenger: &mut Ch,
    ) -> (F, BasefoldProof<F, M>)
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>
            + Clone
            + Send
            + Sync,
    {
        assert_eq!(etas.len(), pd.polys.len(), "one eta per batched poly");
        let combine = |vs: &[Vec<F>]| -> Vec<F> {
            let mut out = F::zero_vec(vs[0].len());
            for (v, &eta) in vs.iter().zip(etas) {
                out.par_iter_mut()
                    .zip(v.par_iter())
                    .for_each(|(o, &x)| *o += eta * x);
            }
            out
        };

        let (v, sumcheck, fold_roots, fold_data) =
            self.prove_folding(combine(&pd.polys), combine(&pd.codewords), z, challenger);
        let pow_witness = tracing::info_span!("basefold_grind").in_scope(|| self.grind(challenger));
        let indices = self.query_indices(challenger);
        let queries = tracing::info_span!("basefold_query_open")
            .in_scope(|| self.collect_queries(&indices, &pd.data, &fold_data));
        let claimed = sumcheck[0][0] + sumcheck[0][1];
        (claimed, BasefoldProof { sumcheck, fold_roots, final_value: v, pow_witness, queries })
    }

    /// The fused eq-sumcheck / codeword-folding engine shared by both opening
    /// paths.
    ///
    /// Drives the transcript: per round observe `h = [h(0), h(1), h(2)]` and
    /// sample `r`; commit and observe each folded oracle `O_1..O_{d-1}`; then
    /// observe the final value. The caller adds the query phase.
    fn prove_folding<Ch>(
        &self,
        mut f: Vec<F>,
        mut cw: Vec<F>,
        z: &[F],
        challenger: &mut Ch,
    ) -> FoldingOutput<F, M>
    where
        Ch: FieldChallenger<F> + CanObserve<M::Commitment>,
    {
        let d = self.num_vars();
        assert_eq!(z.len(), d);
        assert_eq!(f.len(), 1 << d);
        assert_eq!(cw.len(), self.code.codeword_len());

        let mut e = eq_table(z);
        let mut sumcheck = Vec::with_capacity(d);
        let mut fold_roots = Vec::with_capacity(d.saturating_sub(1));
        let mut fold_data: Vec<M::ProverData<RowMajorMatrix<F>>> = Vec::new();

        for i in 0..d {
            let half = f.len() / 2;
            let (fl, fr) = f.split_at(half);
            let (el, er) = e.split_at(half);
            // One parallel pass computes all three round evaluations; three
            // separate reductions would read the same `O(N)` data three times.
            let h = el
                .par_iter()
                .zip(er.par_iter())
                .zip(fl.par_iter().zip(fr.par_iter()))
                .map(|((&a0, &a1), (&b0, &b1))| {
                    Sums([a0 * b0, a1 * b1, (F::TWO * a1 - a0) * (F::TWO * b1 - b0)])
                })
                .sum::<Sums<F, 3>>()
                .0;
            challenger.observe_algebra_slice(&h);
            let r: F = challenger.sample_algebra_element();
            sumcheck.push(h);

            f = bind_msb(&f, r);
            e = bind_msb(&e, r);
            cw = tracing::info_span!("basefold_fold").in_scope(|| self.code.fold(&cw, d - i, r));
            if i + 1 < d {
                let (root, data) = tracing::info_span!("basefold_fold_merkle")
                    .in_scope(|| self.mmcs.commit_matrix(pair_matrix(&cw)));
                challenger.observe(root.clone());
                fold_roots.push(root);
                fold_data.push(data);
            }
        }

        let v = f[0];
        challenger.observe_algebra_element(v);
        (v, sumcheck, fold_roots, fold_data)
    }

    /// Open the fold-pair chain at each sampled index, one **pruned** multi-opening
    /// per round: round `0` from the committed oracle, rounds `1..d` from the
    /// folded oracles.
    fn collect_queries(
        &self,
        indices: &[usize],
        o0_data: &M::ProverData<RowMajorMatrix<F>>,
        fold_data: &[M::ProverData<RowMajorMatrix<F>>],
    ) -> Vec<RoundOpening<F, M>> {
        let d = self.num_vars();
        let n0 = self.code.codeword_len();
        (0..d)
            .map(|i| {
                let n_ip1 = n0 >> (i + 1); // rows of O_i's pair-matrix
                let rows: Vec<usize> = indices.iter().map(|&s| s % n_ip1).collect();
                let data = if i == 0 { o0_data } else { &fold_data[i - 1] };
                let (opened_values, proof) = self.mmcs.open_multi_batch(&rows, data);
                RoundOpening { opened_values, proof }
            })
            .collect()
    }

    /// Verify `f(z) = y` against the commitment to `O_0` (standalone: absorbs the commitment).
    pub fn verify<Ch>(
        &self,
        commitment: &M::Commitment,
        z: &[F],
        y: F,
        proof: &BasefoldProof<F, M>,
        challenger: &mut Ch,
    ) -> Result<(), &'static str>
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>,
    {
        self.verify_with(true, commitment, z, y, proof, challenger)
    }

    /// Verify `g(z) = y` for `g = Σ_j etas[j]·f_j` against a batch commitment
    /// made by [`Self::commit_many`].
    pub fn verify_batch<Ch>(
        &self,
        commitment: &M::Commitment,
        etas: &[F],
        z: &[F],
        y: F,
        proof: &BasefoldProof<F, M>,
        challenger: &mut Ch,
    ) -> Result<(), &'static str>
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>,
    {
        self.verify_core(true, commitment, etas, z, y, proof, challenger)
    }

    /// Batched counterpart of [`Self::verify_opening`]: recovers the value the
    /// opening binds `z` to, without re-absorbing the commitment (an outer
    /// protocol already did).
    pub fn verify_batch_opening<Ch>(
        &self,
        commitment: &M::Commitment,
        etas: &[F],
        z: &[F],
        proof: &BasefoldProof<F, M>,
        challenger: &mut Ch,
    ) -> Result<F, &'static str>
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>,
    {
        if proof.sumcheck.is_empty() {
            return Err("empty sumcheck");
        }
        let v = proof.sumcheck[0][0] + proof.sumcheck[0][1];
        self.verify_core(false, commitment, etas, z, v, proof, challenger)?;
        Ok(v)
    }

    /// Single-polynomial verify. `observe_commitment = false` for the adapter path, where the outer
    /// multi-stark verifier has already absorbed the commitment (re-absorbing would desync
    /// Fiat–Shamir).
    fn verify_with<Ch>(
        &self,
        observe_commitment: bool,
        commitment: &M::Commitment,
        z: &[F],
        y: F,
        proof: &BasefoldProof<F, M>,
        challenger: &mut Ch,
    ) -> Result<(), &'static str>
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>,
    {
        self.verify_core(observe_commitment, commitment, &[F::ONE], z, y, proof, challenger)
    }

    /// The verify engine. `etas` describes the round-0 oracle: a single
    /// `[ONE]` for a plain commitment, or one weight per batched polynomial —
    /// in which case the round-0 folding pair is the `etas`-combination of the
    /// pairs opened from the batch tree. Rounds `1..d` are identical either way,
    /// since the prover folded the combined codeword.
    // One argument over the lint's threshold; every one is load-bearing and the
    // single-poly path funnels through here too.
    #[allow(clippy::too_many_arguments)]
    fn verify_core<Ch>(
        &self,
        observe_commitment: bool,
        commitment: &M::Commitment,
        etas: &[F],
        z: &[F],
        y: F,
        proof: &BasefoldProof<F, M>,
        challenger: &mut Ch,
    ) -> Result<(), &'static str>
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>,
    {
        let d = self.num_vars();
        if etas.is_empty() {
            return Err("empty batch");
        }
        if proof.sumcheck.len() != d || proof.fold_roots.len() != d.saturating_sub(1) {
            return Err("wrong proof shape");
        }
        if observe_commitment {
            challenger.observe(commitment.clone());
        }

        // Replay the sumcheck + fold-commit transcript, deriving r_i.
        let mut claim = y;
        let mut rs = Vec::with_capacity(d);
        for i in 0..d {
            let h = proof.sumcheck[i];
            if h[0] + h[1] != claim {
                return Err("sumcheck round check failed");
            }
            challenger.observe_algebra_slice(&h);
            let r: F = challenger.sample_algebra_element();
            claim = interp_deg2(h, r);
            rs.push(r);
            if i + 1 < d {
                challenger.observe(proof.fold_roots[i].clone());
            }
        }
        if claim != eval_eq(z, &rs) * proof.final_value {
            return Err("final sumcheck ≠ eq(z,r)·v");
        }
        challenger.observe_algebra_element(proof.final_value);

        // COST ABLATION ONLY — never enable in a real verifier. Stops after the
        // sumcheck replay, so the guest's instruction count isolates that phase.
        #[cfg(feature = "ablate-after-sumcheck")]
        return Ok(());

        #[allow(unreachable_code)]
        {
        self.check_grind(proof.pow_witness, challenger)?;
        let indices = self.query_indices(challenger);
        if proof.queries.len() != d {
            return Err("wrong per-round opening count");
        }
        let n0 = self.code.codeword_len();
        let od = self.code.base_encode(proof.final_value); // O_d = v·g0

        // One pruned multi-opening per round covers every query, so the Merkle
        // check is per round; the fold chain is then walked per query.
        for (i, round) in proof.queries.iter().enumerate() {
            let n_ip1 = n0 >> (i + 1);
            let width = if i == 0 { 2 * etas.len() } else { 2 };
            let rows: Vec<usize> = indices.iter().map(|&s| s % n_ip1).collect();
            if round.opened_values.len() != rows.len() {
                return Err("wrong opened-row count");
            }
            if round
                .opened_values
                .iter()
                .any(|v| v.len() != 1 || v[0].len() != width)
            {
                return Err("opened row has the wrong width");
            }
            let dims = [Dimensions { width, height: n_ip1 }];
            self.mmcs
                .verify_multi_batch(
                    if i == 0 { commitment } else { &proof.fold_roots[i - 1] },
                    &dims,
                    &rows,
                    &round.opened_values,
                    &round.proof,
                )
                .map_err(|_| "merkle verify failed")?;
        }

        // COST ABLATION ONLY — stops after Merkle verification, before the
        // per-query fold-consistency walk.
        #[cfg(feature = "ablate-after-merkle")]
        return Ok(());

        #[allow(unreachable_code)]
        for (qi, &s) in indices.iter().enumerate() {
            let mut prev_folded: Option<F> = None;
            for i in 0..d {
                let n_ip1 = n0 >> (i + 1);
                let row = s % n_ip1;
                let vals = &proof.queries[i].opened_values[qi][0];
                // Round 0 reconstructs the combined oracle's pair from the
                // batched pairs; later rounds already are the combined oracle.
                let (a, b) = if i == 0 {
                    etas.iter().enumerate().fold((F::ZERO, F::ZERO), |(a, b), (j, &eta)| {
                        (a + eta * vals[2 * j], b + eta * vals[2 * j + 1])
                    })
                } else {
                    (vals[0], vals[1])
                };

                if let Some(pf) = prev_folded {
                    // The value produced by the previous fold sits at position (s mod N_i) in O_i,
                    // which is either the low (a) or high (b) sibling of this opened pair.
                    let n_i = n0 >> i;
                    let pos = s % n_i;
                    let expected = if pos < n_ip1 { a } else { b };
                    if pf != expected {
                        return Err("fold-chain consistency failed");
                    }
                }
                prev_folded = Some(self.code.fold_pair(a, b, d - i, row, rs[i]));
            }
            // Last fold must land on O_d = v·g0.
            let last = prev_folded.expect("d >= 1");
            let last_row = s % od.len();
            if last != od[last_row] {
                return Err("final oracle check failed");
            }
        }
        Ok(())
        }
    }

    /// Integration-facing verify (mirrors `TensorPcs::verify_opening`): does NOT take an external
    /// claimed value and does NOT re-absorb the commitment. Recovers and returns `v = P(z)` — the
    /// value the opening binds `z` to — for the outer adapter to check against its batched claim.
    /// `v` is read as `h_0(0)+h_0(1)` and bound to the committed `P` by the sumcheck chain + fold.
    pub fn verify_opening<Ch>(
        &self,
        commitment: &M::Commitment,
        z: &[F],
        proof: &BasefoldProof<F, M>,
        challenger: &mut Ch,
    ) -> Result<F, &'static str>
    where
        Ch: FieldChallenger<F>
            + CanObserve<M::Commitment>
            + CanSampleBits<usize>
            + GrindingChallenger<Witness = F>,
    {
        if proof.sumcheck.is_empty() {
            return Err("empty sumcheck");
        }
        let v = proof.sumcheck[0][0] + proof.sumcheck[0][1];
        // Adapter path: do NOT re-absorb the commitment (the outer verifier already did).
        self.verify_with(false, commitment, z, v, proof, challenger)?;
        Ok(v)
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use p3_challenger::DuplexChallengerPF;
    use p3_ed25519::{
        ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519, Poseidon2Ed25519,
    };
    use p3_field::PrimeCharacteristicRing;
    use p3_merkle_tree::MerkleTreeMmcs;
    use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    use super::*;
    use crate::pcs::eval_multilinear;

    // Match the local eddsa-example config: Poseidon2 over Ed25519 (width 3), DuplexChallengerPF.
    type F = Ed25519;
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    const OUT: usize = 1;
    type Perm = Poseidon2Ed25519<WIDTH>;
    type MyHash = PaddingFreeSponge<Perm, WIDTH, RATE, OUT>;
    type MyCompress = TruncatedPermutation<Perm, 2, OUT, WIDTH>;
    type MyMmcs = MerkleTreeMmcs<Ed25519, Ed25519, MyHash, MyCompress, 2, OUT>;
    type Ch = DuplexChallengerPF<Ed25519, Perm, WIDTH, RATE>;

    fn perm() -> Perm {
        // Deterministic (ChaCha) so every construction is identical.
        let mut rng = ChaCha8Rng::seed_from_u64(0xED25519_C0FFEE);
        Poseidon2Ed25519::<WIDTH>::new_from_rng(
            ED25519_POSEIDON2_FULL_ROUNDS,
            ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
            &mut rng,
        )
    }
    fn mmcs() -> MyMmcs {
        let p = perm();
        MyMmcs::new(MyHash::new(p.clone()), MyCompress::new(p), 0)
    }
    fn challenger() -> Ch {
        Ch::new(perm())
    }

    fn pcs(d: usize, q: usize) -> BasefoldPcs<F, MyMmcs> {
        BasefoldPcs::new(RandomFoldableCode::<F>::new(d, 4, 0xBA5E_F01D), mmcs(), q)
    }

    #[test]
    fn honest_open_verifies() {
        for d in 2..=8 {
            let p = pcs(d, 12);
            let mut rng = ChaCha8Rng::seed_from_u64(d as u64);
            let evals: Vec<F> = (0..(1 << d)).map(|_| rng.random::<F>()).collect();
            let z: Vec<F> = (0..d).map(|_| rng.random::<F>()).collect();
            let y = eval_multilinear(&evals, &z);

            let mut ch = challenger(); // commit + open share one transcript
            let (commitment, pd) = p.commit(evals, &mut ch);
            let (_v, proof) = p.open(&pd, &z, &mut ch);
            p.verify(&commitment, &z, y, &proof, &mut challenger())
                .unwrap_or_else(|e| panic!("honest reject at d={d}: {e}"));
        }
    }

    #[test]
    fn wrong_value_rejected() {
        let p = pcs(6, 12);
        let mut rng = ChaCha8Rng::seed_from_u64(1);
        let evals: Vec<F> = (0..64).map(|_| rng.random::<F>()).collect();
        let z: Vec<F> = (0..6).map(|_| rng.random::<F>()).collect();
        let y = eval_multilinear(&evals, &z);
        let mut ch = challenger();
        let (commitment, pd) = p.commit(evals, &mut ch);
        let (_v, proof) = p.open(&pd, &z, &mut ch);
        assert!(p.verify(&commitment, &z, y + F::ONE, &proof, &mut challenger()).is_err());
    }

    #[test]
    fn tampered_query_rejected() {
        let p = pcs(6, 20);
        let mut rng = ChaCha8Rng::seed_from_u64(2);
        let evals: Vec<F> = (0..64).map(|_| rng.random::<F>()).collect();
        let z: Vec<F> = (0..6).map(|_| rng.random::<F>()).collect();
        let y = eval_multilinear(&evals, &z);
        let mut ch = challenger();
        let (commitment, pd) = p.commit(evals, &mut ch);
        let (_v, mut proof) = p.open(&pd, &z, &mut ch);
        // Corrupt a claimed final value -> final-oracle + eq(z,r)·v checks must fail.
        proof.final_value += F::ONE;
        assert!(p.verify(&commitment, &z, y, &proof, &mut challenger()).is_err());
    }

    fn rand_batch(rng: &mut ChaCha8Rng, k: usize, d: usize) -> Vec<Vec<F>> {
        (0..k)
            .map(|_| (0..(1 << d)).map(|_| rng.random::<F>()).collect())
            .collect()
    }

    /// A batched opening proves the eta-combination of every committed poly, at
    /// one point, with a single proof.
    #[test]
    fn batch_open_verifies() {
        for (k, d) in [(1, 4), (2, 5), (4, 6), (5, 3)] {
            let p = pcs(d, 12);
            let mut rng = ChaCha8Rng::seed_from_u64(100 + k as u64);
            let polys = rand_batch(&mut rng, k, d);
            let etas: Vec<F> = (0..k).map(|_| rng.random::<F>()).collect();
            let z: Vec<F> = (0..d).map(|_| rng.random::<F>()).collect();

            // The claim is the eta-combination of the individual evaluations.
            let y: F = polys
                .iter()
                .zip(&etas)
                .map(|(f, &eta)| eta * eval_multilinear(f, &z))
                .sum();

            let mut ch = challenger();
            let (commitment, pd) = p.commit_many(polys, &mut ch);
            let (v, proof) = p.open_batch(&pd, &etas, &z, &mut ch);
            assert_eq!(v, y, "batched opening value wrong at k={k} d={d}");

            p.verify_batch(&commitment, &etas, &z, y, &proof, &mut challenger())
                .unwrap_or_else(|e| panic!("honest batch reject at k={k} d={d}: {e}"));
        }
    }

    /// `k = 1` must agree with the single-polynomial path bit for bit: same
    /// commitment, same proof shape, same accepted value.
    #[test]
    fn batch_of_one_matches_single() {
        let p = pcs(6, 12);
        let mut rng = ChaCha8Rng::seed_from_u64(7);
        let evals: Vec<F> = (0..64).map(|_| rng.random::<F>()).collect();
        let z: Vec<F> = (0..6).map(|_| rng.random::<F>()).collect();

        let mut ch_single = challenger();
        let (c_single, pd_single) = p.commit(evals.clone(), &mut ch_single);
        let (v_single, _) = p.open(&pd_single, &z, &mut ch_single);

        let mut ch_batch = challenger();
        let (c_batch, pd_batch) = p.commit_many(vec![evals.clone()], &mut ch_batch);
        let (y_batch, proof_batch) = p.open_batch(&pd_batch, &[F::ONE], &z, &mut ch_batch);

        assert_eq!(c_single, c_batch, "k=1 commitment must match the single path");
        // `open` returns f(r); `open_batch` returns the claimed value f(z). Both
        // proofs carry the same final folded value.
        assert_eq!(v_single, proof_batch.final_value);
        assert_eq!(y_batch, eval_multilinear(&evals, &z));
    }

    #[test]
    fn batch_rejects_wrong_value_and_wrong_etas() {
        let (k, d) = (3, 6);
        let p = pcs(d, 16);
        let mut rng = ChaCha8Rng::seed_from_u64(8);
        let polys = rand_batch(&mut rng, k, d);
        let etas: Vec<F> = (0..k).map(|_| rng.random::<F>()).collect();
        let z: Vec<F> = (0..d).map(|_| rng.random::<F>()).collect();

        let mut ch = challenger();
        let (commitment, pd) = p.commit_many(polys, &mut ch);
        let (y, proof) = p.open_batch(&pd, &etas, &z, &mut ch);

        assert!(
            p.verify_batch(&commitment, &etas, &z, y + F::ONE, &proof, &mut challenger())
                .is_err(),
            "a wrong claimed value must be rejected"
        );

        // Verifying with different weights checks a different virtual poly, so
        // the round-0 pair no longer folds into the committed chain.
        let mut other = etas.clone();
        other[1] += F::ONE;
        assert!(
            p.verify_batch(&commitment, &other, &z, y, &proof, &mut challenger())
                .is_err(),
            "mismatched etas must be rejected"
        );
    }

    #[test]
    fn batch_rejects_tampered_leaf() {
        let (k, d) = (2, 6);
        let p = pcs(d, 20);
        let mut rng = ChaCha8Rng::seed_from_u64(9);
        let polys = rand_batch(&mut rng, k, d);
        let etas: Vec<F> = (0..k).map(|_| rng.random::<F>()).collect();
        let z: Vec<F> = (0..d).map(|_| rng.random::<F>()).collect();

        let mut ch = challenger();
        let (commitment, pd) = p.commit_many(polys, &mut ch);
        let (y, mut proof) = p.open_batch(&pd, &etas, &z, &mut ch);

        // Corrupt one batched round-0 leaf value; its Merkle path breaks.
        proof.queries[0].opened_values[0][0][0] += F::ONE;
        assert!(
            p.verify_batch(&commitment, &etas, &z, y, &proof, &mut challenger())
                .is_err(),
            "a tampered batched leaf must be rejected"
        );
    }

    /// Grinding must actually gate the query phase: a proof carrying the wrong
    /// witness is rejected, and one carrying the right witness verifies.
    #[test]
    fn grinding_gates_the_query_phase() {
        const POW: usize = 8;
        let p = BasefoldPcs::new_with_grinding(
            RandomFoldableCode::<F>::new(6, 4, 0xBA5E_F01D),
            mmcs(),
            12,
            POW,
        );
        let mut rng = ChaCha8Rng::seed_from_u64(21);
        let evals: Vec<F> = (0..64).map(|_| rng.random::<F>()).collect();
        let z: Vec<F> = (0..6).map(|_| rng.random::<F>()).collect();
        let y = eval_multilinear(&evals, &z);

        let mut ch = challenger();
        let (commitment, pd) = p.commit(evals, &mut ch);
        let (_v, mut proof) = p.open(&pd, &z, &mut ch);
        p.verify(&commitment, &z, y, &proof, &mut challenger())
            .expect("a ground proof must verify");

        proof.pow_witness += F::ONE;
        assert!(
            p.verify(&commitment, &z, y, &proof, &mut challenger()).is_err(),
            "a bad proof-of-work witness must be rejected"
        );
    }
}
