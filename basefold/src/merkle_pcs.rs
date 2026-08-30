//! M2-auth: the **succinct** Basefold PCS — Merkle-committed folded oracles + a `q_bf`-query phase
//! driven by a real Fiat–Shamir challenger. This is the polylog-verifier version (the M2-core in
//! `pcs.rs` sent full oracles; here the verifier only checks `q_bf` query paths).
//!
//! Single field (`Challenge = Val`), matching alp's eddsa config (`EF = Ed25519 = F`).

use alloc::vec::Vec;

use p3_challenger::{CanObserve, CanSampleBits, FieldChallenger};
use p3_commit::{BatchOpening, BatchOpeningRef, Mmcs};
use p3_field::Field;
use p3_matrix::Dimensions;
use p3_matrix::dense::RowMajorMatrix;
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
    /// Per query `q`, per round `i` (0..d): the opened folding-pair of `O_i` + Merkle proof.
    pub queries: Vec<Vec<BatchOpening<F, M>>>,
}

/// Reshape a codeword `cw` (length `n`) into an `(n/2) × 2` matrix whose row `j` is the folding
/// pair `[cw[j], cw[j+n/2]]`, so one Merkle path opens both siblings of a fold.
fn pair_matrix<F: Field>(cw: &[F]) -> RowMajorMatrix<F> {
    let half = cw.len() / 2;
    let mut flat = Vec::with_capacity(cw.len());
    for j in 0..half {
        flat.push(cw[j]);
        flat.push(cw[j + half]);
    }
    RowMajorMatrix::new(flat, 2)
}

/// The succinct Basefold PCS over a foldable code `code` and a Merkle commitment scheme `mmcs`.
pub struct BasefoldPcs<F: Field, M: Mmcs<F>> {
    pub code: RandomFoldableCode<F>,
    pub mmcs: M,
    pub num_queries: usize,
}

impl<F: Field, M: Mmcs<F>> BasefoldPcs<F, M> {
    pub fn new(code: RandomFoldableCode<F>, mmcs: M, num_queries: usize) -> Self {
        Self { code, mmcs, num_queries }
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

    /// Prove `f(z) = y`.
    pub fn open<Ch>(&self, pd: &BasefoldProverData<F, M>, z: &[F], challenger: &mut Ch) -> (F, BasefoldProof<F, M>)
    where
        Ch: FieldChallenger<F> + CanObserve<M::Commitment> + CanSampleBits<usize>,
    {
        let d = self.num_vars();
        assert_eq!(z.len(), d);

        let mut f = pd.evals.clone();
        let mut e = eq_table(z);
        let mut cw = pd.o0_cw.clone();

        let mut sumcheck = Vec::with_capacity(d);
        let mut fold_roots = Vec::with_capacity(d.saturating_sub(1));
        let mut fold_data: Vec<M::ProverData<RowMajorMatrix<F>>> = Vec::new();

        for i in 0..d {
            let half = f.len() / 2;
            let (fl, fr) = f.split_at(half);
            let (el, er) = e.split_at(half);
            let h0: F = el.iter().zip(fl).map(|(&a, &b)| a * b).sum();
            let h1: F = er.iter().zip(fr).map(|(&a, &b)| a * b).sum();
            let h2: F = (0..half).map(|x| (F::TWO * er[x] - el[x]) * (F::TWO * fr[x] - fl[x])).sum();
            let h = [h0, h1, h2];
            challenger.observe_algebra_slice(&h);
            let r: F = challenger.sample_algebra_element();
            sumcheck.push(h);

            f = bind_msb(&f, r);
            e = bind_msb(&e, r);
            cw = self.code.fold(&cw, d - i, r);
            if i + 1 < d {
                let (root, data) = self.mmcs.commit_matrix(pair_matrix(&cw));
                challenger.observe(root.clone());
                fold_roots.push(root);
                fold_data.push(data);
            }
        }

        let v = f[0];
        challenger.observe_algebra_element(v);

        let indices = self.query_indices(challenger);
        let n0 = self.code.codeword_len();
        let queries = indices
            .iter()
            .map(|&s| {
                (0..d)
                    .map(|i| {
                        let n_ip1 = n0 >> (i + 1); // rows of O_i's pair-matrix
                        let row = s % n_ip1;
                        if i == 0 {
                            self.mmcs.open_batch(row, &pd.o0_data)
                        } else {
                            self.mmcs.open_batch(row, &fold_data[i - 1])
                        }
                    })
                    .collect()
            })
            .collect();

        (v, BasefoldProof { sumcheck, fold_roots, final_value: v, queries })
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
        Ch: FieldChallenger<F> + CanObserve<M::Commitment> + CanSampleBits<usize>,
    {
        self.verify_with(true, commitment, z, y, proof, challenger)
    }

    /// Core verify. `observe_commitment = false` for the adapter path, where the outer multi-stark
    /// verifier has already absorbed the commitment (re-absorbing would desync Fiat–Shamir).
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
        Ch: FieldChallenger<F> + CanObserve<M::Commitment> + CanSampleBits<usize>,
    {
        let d = self.num_vars();
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

        let indices = self.query_indices(challenger);
        if proof.queries.len() != indices.len() {
            return Err("wrong query count");
        }
        let n0 = self.code.codeword_len();
        let od = self.code.base_encode(proof.final_value); // O_d = v·g0

        for (qi, &s) in indices.iter().enumerate() {
            let opening = &proof.queries[qi];
            if opening.len() != d {
                return Err("wrong per-query round count");
            }
            let mut prev_folded: Option<F> = None;
            for i in 0..d {
                let n_ip1 = n0 >> (i + 1);
                let row = s % n_ip1;
                let root = if i == 0 { commitment } else { &proof.fold_roots[i - 1] };
                let dims = [Dimensions { width: 2, height: n_ip1 }];
                self.mmcs
                    .verify_batch(root, &dims, row, BatchOpeningRef::from(&opening[i]))
                    .map_err(|_| "merkle verify failed")?;
                let vals = &opening[i].opened_values[0];
                if vals.len() != 2 {
                    return Err("opened row not a pair");
                }
                let (a, b) = (vals[0], vals[1]);

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
        Ch: FieldChallenger<F> + CanObserve<M::Commitment> + CanSampleBits<usize>,
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
}
