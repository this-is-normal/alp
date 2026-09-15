# Quasar PCS — implementation plan

**Goal:** add Quasar as a third multilinear PCS in this fork, next to `brakedown/` (tensor) and
`basefold/` (foldable-code). Quasar = **QA-code commit front end** (Brakedown-shaped, but with a
WHT encoder) + **code switching into BaseFold**, giving a linear-ish prover with a *polylog*
verifier and proof.

Paper: **"Quasar: A Field-Agnostic Polynomial Commitment Scheme with Polylogarithmic Verification
from Quasi-Abelian Codes"** — Jia, Li, Xing, Yao, Yuan (eprint 2026/1839).
Local PDF: `../2026-1839.pdf`; extracted text: scratchpad `quasar.txt`.

**Authors' code is available and mirrored locally** (see §1) — the QA encoder, the security-parameter
formulas and the two-layer sumcheck are all portable, so most of this plan is *porting*, not
inventing.

---

## 0. Why this scheme belongs in *this* fork

Quasar is field-agnostic and FFT-free — the same reason `brakedown/` and `basefold/` exist here.
Better: the QA distance bound *improves* with field size, so the ed25519 base field (255 bits) is
the best case, not a handicap. Ported bound (`qabase_distance_lower_bound`, λ_dist = 103, c = 2):

| field bits | δ_QA @ N₂=2¹² | δ_QA @ N₂=2²⁰ | t (columns opened) |
|---|---|---|---|
| 31 (BabyBear/KoalaBear) | **0** — unusable | 0 | — |
| 64 (Goldilocks) | 0.250 | 0.230 | 821–896 |
| 127 (paper's Mersenne127) | 0.445 | 0.414 | 445–482 |
| **255 (ed25519)** | **0.473** | **0.457** | **417–433** |

At inverse rate c = 4 the bound jumps to δ ≈ 0.73 (t ≈ 256). So: **ed25519 is the bring-up field**,
and 31-bit fields are out of scope for QA by construction.

---

## 1. Sources and what we reuse

### Authors' artifact (mirrored)

`/home/ubuntu/code/alp_agents/quasar-ref/` — full mirror of
`https://anonymous.4open.science/r/Quasar-33F6` (124 files, 1.9 MB; downloaded 2026-09-06).
It is a fork of `hadasz/plonkish_basefold` (which already contains BaseFold + Brakedown).

| artifact file | what it gives us |
|---|---|
| `plonkish_backend/src/pcs/multilinear/quasar.rs` (2455 L) | **the reference Quasar**: security config, WHT, QA encode, selector sumcheck (`prove_merged_column_consistency`), two-layer GKR (`prove_qabase_two_layer_gkr_batch`), full open/verify, 5 tests |
| `.../qabase.rs` (5308 L) | earlier/companion superset (aux commitments, witness-carrying encode) — useful cross-reference |
| `.../qapcs.rs` (1273 L) | the Brakedown-style QA baseline (**QAPCS**), incl. matrix-shape scoring |
| `.../basefold.rs`, `brakingbase.rs`, `brakedown.rs` | the baselines they measured against |
| `benchmark/benches/quasar_bench.rs`, `matmul_bench.rs` | CLI/param plumbing, the exact experiment shapes |
| `qa_gpu_overlay/cuda/qa_mersenne127.cu` (32 KB) | CUDA WHT encoder + GPU column hashing (later milestone) |

No `bench_data/` CSVs in the artifact — the only reference numbers are the paper's Table 4 / Fig. 2 /
Table 3.

### This repo

| we reuse | for |
|---|---|
| `brakedown/`: `LinearCode`/`LinearCodeExt`, `TensorPcs`, `BrakedownMlPcs`, `ReedSolomonCode` | the QA code plugs straight into these traits → **QAPCS for free** |
| `basefold/`: `RandomFoldableCode`, `BasefoldPcs` (Merkle + query phase), `BasefoldMlPcs` | the code-switching target (black-box BaseFold in the paper) |
| `sumcheck/`: `ProductPolynomial`, `generic_degree`, `Constraint`/weights, `PrescribedPointPcs` | the selector sumcheck + the two QA-encoding sumchecks |
| `security/` (`p3-security`) | the designated home for the QA distance/query math |
| `merkle-tree`, `challenger` (`DuplexChallengerCF`), `symmetric` (`SplitFieldHasher`) | ed25519 data + BabyBear hash stack, already the house config |
| `pcs-profile/` (`ProfileablePcs`, `PcsMetrics`) + `stark-stats` | **the estimator-validation harness** (§5) |
| `commit/`: `MultilinearPcs` | the trait Quasar must implement |

### What is genuinely new code

The code-switching layer: selector sumcheck over the sampled columns, the two-layer WHT/pointwise
sumcheck, and **batched** BaseFold opening (our `BasefoldPcs` opens one vector at one point today).

---

## 2. The protocol, as it will be implemented

Shape: `N = 2ⁿ = N₁·N₂`, `N₁ = 2ⁿ¹` rows (small, 2⁶–2⁸), `N₂` columns, QA inverse rate `c`.

**QA code (rank-1 systematic, index c).** Messages are elements of `R = F[X₁..X_{n₂}]/(Xᵢ²−1)`;
multiplication in `R` = pointwise multiply after Walsh–Hadamard transform. Encoding a row `v`:

```
v'   = WHT(v)                       # H = [[1,1],[1,-1]], H^{⊗n₂}
u_i' = (2^{-n₂} a_i') ⋆ v'          # public scaled vectors, a_1' = 1
u_i  = WHT(u_i')                    # codeword = (u_1 ‖ … ‖ u_c), u_1 = v
```

Cost per row: `2c·N₂ log N₂` adds, `cN₂` muls — **no field muls in the WHT itself**, which is why
commit is fast (paper Table 1: `N` muls, `2N log N` adds, `4N−1` hashes at rate ½).

**Commit** = reshape to `M ∈ F^{N₁×N₂}`, QA-encode each row → `C ∈ F^{N₁×cN₂}`, Merkle-commit the
*columns*. Identical in shape to `TensorPcs::commit`.

**Eval** (both branches run in parallel — proximity with random `r`, evaluation with `z_L`):

1. `v = rᵀM`, `p = z_LᵀM`; `u = QA.Enc(v)`, `q = QA.Enc(p)`.
2. BaseFold-commit the `c` blocks of `u` and of `q` (length `N₂` each).
3. Sample `t = Θ(λ)` columns of `C`; open them with Merkle paths.
4. **Selector sumcheck** (`log cN₂` rounds, degree 2): with sparse `h̃(⟨i_j⟩) = τ^j`,
   `Σ_x h̃(x)ũ(x) = Σ_j τ^j⟨r, C(·,i_j)⟩` and the same for `q̃` / `z_L`. Terminal: block evals at `e₁`.
5. **QA-encoding check, layer 1** (`n₂` rounds, degree 3):
   `Σ_i α^{i-2} ũ_i(γ) = Σ_x H̃(γ,x)·(Σ_i α^{i-2} ã_i'(x))·ṽ'(x)`, batched over both branches.
   `H̃(γ,·)` is materialized on the hypercube (a tensor product, `∏_j(1−2γ_j x_j)`); the `ã_i'` are
   *preprocessed* BaseFold commitments; `v' = WHT(v)` is a prover-side witness, **never committed**.
6. **QA-encoding check, layer 2** (`n₂` rounds, degree 2): substitute `ṽ'(r) = Σ_y H̃(r,y)ṽ(y)`,
   reducing to an evaluation of the *committed input block* `v` (= `u₁`). This is what kills the
   commitment to `v'`.
7. **One batched BaseFold opening** for every claim from steps 4/5/6 plus the final `q̃₁(z_R) = y`.
   The verifier evaluates `H̃` and `h̃` itself in `O(log N)`.

Verifier: `O(λ log N)` for the columns + `O(λ log² N₂)` for the batched BaseFold ⇒ `O(λ log² N)`.

---

## 3. Milestones

Each milestone ends with a green test command and, from M2 on, a number recorded in the status log.

### M0 — recon *(done)*
- Artifact mirrored; paper's construction/efficiency/parameter sections read.
- **δ_QA formula ported to Python and validated against the paper's Table 4: all 10 rows match
  exactly** (δ to 6 dp and `t` exactly, k = 20…29). Scratchpad `dist_check.py`.

### M1 — `p3-security`: the QA soundness math
`security/src/qa.rs`, mirroring `quasar.rs:74-153`:
- `qa_gp(δ, field_bits)`, `qa_distance_failure_log2(δ, n₂, c, field_bits)` (bisection target),
  `qa_distance_lower_bound(...)`, `qa_num_queries(δ, λ) = ⌈λ / −log₂(1−δ/3)⌉`,
  and `QaRegime` in the style of the crate's other regimes.
- **Test:** the 10 rows of paper Table 4 (Mersenne127, c=2, λ_dist=103) as a hard-coded fixture,
  plus the ed25519/`c=4` values, plus "31-bit ⇒ no admissible δ" returning an error not a panic.
- `no_std` + `libm` for `log2`/`exp_m1` (the crate already depends on libm).

### M2 — the QA code + QAPCS baseline (the cheap win)
New workspace member **`qa-code/`** (`p3-qa-code`), depending on `p3-brakedown` for the
`LinearCode`/`LinearCodeExt` traits. (If that dependency direction grates later, the two traits lift
out of `brakedown/src/lib.rs` into a `linear-code/` crate and get re-exported from `brakedown` — no
call-site changes anywhere. Not worth doing pre-emptively.)
- `wht<F>(&mut [F])` and a maybe-rayon parallel variant, over `RowMajorMatrix` columns (the
  `LinearCode` convention: `message_len` rows × any width, each column an independent message).
- `QaCode<F>`: `new_random(n₂, c, rng)` storing the scaled `2^{-n₂}a_i'`; `impl LinearCode<F>` +
  `impl LinearCodeExt<F, EF>`.
- **Tests:** systematic (`u₁ == v`); `encode` == naive dyadic-matrix multiply for `n₂ ≤ 6`;
  `F`-linearity; `encode_ext` agrees with `encode` on embedded base values; `c ∈ {2,4}`.
- **Free deliverable:** `TensorPcs<QaCode>` / `BrakedownMlPcs<…QaCode…>` = **QAPCS**. Reuse
  `brakedown/tests/tensor_pcs_ed25519.rs` and `examples/pcs_cost_sweep_ed25519.rs` verbatim with the
  code swapped, with `t` from M1.
- **Numbers to record:** QA vs Brakedown vs RS encode time, and QAPCS commit/open/verify at
  2¹⁶–2²², over ed25519 on vultr3.

### M3 — `quasar/` crate: params + commit phase
- `QuasarParams`: `n`, `n₁`, `c`, `λ`, derived `δ_QA`, `t`; a shape chooser (port the scoring proxy
  from `qapcs.rs:297-366`, extended to Quasar's polylog proof so it picks small `N₁`).
- `QuasarPcs::commit` — reshape, QA row-encode, column-Merkle via `Mmcs`. Reuse `TensorPcs`'s
  commit path rather than re-writing it; keep `M` and `C` in prover data.
- Setup: `RandomFoldableCode` for the length-`N₂` BaseFold instances + **preprocessed BaseFold
  commitments to the `a_i'`**.
- **Test:** commit is deterministic and equals the QAPCS root for the same code/shape.

### M4 — selector sumcheck (sampled-column consistency)
- Sample `t` column indices from the challenger (`CanSampleBits`), open them from the MMCS.
- Sparse `h̃` (length `cN₂`, `t` nonzeros) and the degree-2 sumcheck against the concatenated
  `ũ`/`q̃`; batch the two branches with one challenge.
- Verifier side: `Σ_j τ^j⟨r, C(·,i_j)⟩`, path checks, and `h̃(e₁)` by direct sparse evaluation.
- **Tests:** honest transcript accepts; tampered opened column rejects; wrong `τ`-weighting rejects.

### M5 — the two-layer QA-encoding sumcheck
- Layer 1: degree-3 product `v' · H̃(γ,·) · E_i` batched over blocks × branches
  (ref: `quasar.rs:1353-1560`). Use `sumcheck::generic_degree` (or `ProductPolynomial` with the
  `H̃` factor carried as a weight, which avoids materializing it — decide by benchmark).
- Layer 2: degree-2 `H̃(outer,·) · (Σ_b w_b v_b)`, terminal claims on the committed input blocks.
- `H̃` helpers: `hadamard_tensor_mle_eval(a,b) = ∏(1−2a_jb_j)` and its hypercube expansion.
- **Tests:** a standalone "QA-encoding relation" test (prove that a *correct* codeword's blocks
  satisfy the relation, and that a single corrupted output entry is rejected); `c ∈ {2,4}`.

### M6 — batched BaseFold opening (the piece our `basefold/` crate lacks)
All length-`N₂` vectors (`u₁..u_c`, `q₁..q_c`, and the preprocessed `a_i'`) must be opened at several
points with **one** proof, or the proof size blows up (≈ 300 KB per naive opening × ~8).
- Commit the `2c` online blocks as a **multi-column MMCS matrix** so one query index opens all
  blocks at once, then fold the `η`-linear combination leaf-wise.
- Add `BasefoldPcs::{commit_batch, open_batch, verify_batch}`: an `ℓ`-round sumcheck reducing all
  `(point, value)` claims to a common `ζ`, then one BaseFold evaluation of `Σ ηₜ fₜ`
  (paper §2 "batch evaluation protocol", following NST25).
- **Tests:** `open_batch` on `k = 1` matches today's `open`; `k = 5` claims at 3 distinct points
  round-trips; each of (wrong value, tampered leaf, swapped `η`) rejects.

### M7 — end-to-end Quasar
- `QuasarPcs::{open, verify}` wiring M4+M5+M6, plus the final `q̃₁(z_R) = y` check.
- `impl MultilinearPcs` and (for `multi-stark` drop-in next to `BrakedownMlPcs`/`BasefoldMlPcs`)
  `impl PrescribedPointPcs` → `QuasarMlPcs`.
- **Tests:** roundtrip; wrong claimed value; tampered column; tampered BaseFold leaf; corrupted
  encoding witness; shapes `n₁ ∈ {5,6,8}` × `c ∈ {2,4}` × `n ∈ {12,16,20}`.
- Criterion bench + `pcs-profile/src/quasar.rs` (`ProfileablePcs`) so Quasar joins the `sweep` bin.

### M8 — optional / later (PCS-level)
- CUDA WHT encoder on vultr3's RTX 3060, porting `qa_mersenne127.cu` (the paper's biggest lever:
  33–148× on encode; but note commit is ~97% Merkle hashing in our ed25519 config, so the ceiling
  is small until leaf hashing also moves to the GPU).
- `eddsa-example --features quasar` (the real end-to-end signal: does a polylog verifier shrink the
  OpenVM recursion guest?).

---

## 3b. Phase 2 — end-to-end proving

M1–M7 delivered Quasar as a *PCS*. Phase 2 is proving a real AIR with it: `multi-stark` +
`QuasarMlPcs`, ending at `eddsa-example --features quasar` (a full ed25519 signature-verification
STARK), measured against the Brakedown, Ligero and BaseFold backends already in the tree.

### E0 — spike: does the adapter actually run? *(done 2026-09-06)*

`QuasarMlPcs` compiled at M7 but had never been executed. Wiring it into `ed25519-fib-example`
(`config_for_quasar` + `tests/quasar_pcs.rs`) surfaced and fixed **two real gaps**, and both
end-to-end tests now pass:

- `QuasarProverData` and `QuasarMlProverData` needed `Clone` — the multi-stark prover clones prover
  data across opening rounds. (`BasefoldMlProverData` has the same manual impl; ours was missing, so
  the adapter could not have been used by anything.)
- Nothing else. The transcript, the layout sumcheck, the prescribed-point opening and the
  preprocessed-coefficient commitment all line up on the first run.

So the remaining work is integration, parameters and measurement — not protocol debugging.

### E1 — the real target: `eddsa-example --features quasar`

`EdDsaAir` declares a **preprocessed trace**, so the config carries two PCS instances at two
different stacked arities (main width 45, preprocessed width 8) — the one code path the spike did
not exercise.

- `eddsa-example/src/config.rs` currently selects its backend with `#[cfg(not(feature = "basefold"))]`
  style conditions across ~15 sites. A fourth backend makes that unreadable: **restructure into one
  `mod backend_{brakedown,ligero,basefold,quasar}` per PCS with a single re-export**, then add the
  Quasar arm. This is the bulk of E1 and it pays for itself the next time a PCS is added.
- Two Quasar instances mean two QA codes, two BaseFold codes and two preprocessed coefficient
  commitments. Note for later: the *verifier* rebuilds the coefficient commitment at setup
  (an encode + Merkle over `c−1` length-`N₂` vectors); in a deployment it belongs in the verifying
  key as a pinned constant, as the OpenVM recursion work already does for `vk`.
- **Tests:** the `eddsa_brakedown.rs` pair (honest roundtrip + tampered public value) re-run under
  `--features quasar`.

### E2 — parameter selection from the cost model

The spike picks `log_rows = min(3, stacked−1)` by hand. Replace with a real chooser:

- `QuasarParams::choose(stacked_num_variables, field_bits, target_bits)` minimising predicted proof
  bytes (with a switch for verifier cost), mirroring `p3_brakedown::params::optimal_log_rows`.
- Validate it against a brute-force sweep of `quasar_cost_model.py`, which already predicts proof
  bytes to within 1 %. The model says the optimum is shallow and sits at `N₁ = 2⁶–2⁷` at 2²⁵–2²⁹, but
  small traces are a different regime and need their own answer.
- The BaseFold inverse rate (4) and repetition count (288) are currently constants by fiat; derive
  them from a target bit count through `p3-security` the way the column count already is.

### E3 — honest end-to-end security

- At small arities the column count clamps to the codeword and the achieved soundness collapses —
  the fib demo runs at ~8 bits. Add `QuasarParams::achieved_security_bits()` and have every config
  log it at build time (the Ligero path already logs its query count); make the examples state the
  regime they are in rather than implying 100 bits.
- Compose the PCS error with the zerocheck/AIR error via `p3-security::stark::proven_security`, so
  the number quoted is for the *proof*, not for the PCS in isolation.

### E4 — measurement and the crossover question

`ed25519-fib-example/tests/pcs_comparison.rs` (added in E0) prints proof bytes and prove/verify time
per backend for one trace. Extend it to the eddsa AIR and all four backends, and answer the two
questions the early numbers raise:

**PCS-only, ed25519 + BabyBear hash, `N₁ = 2⁶`, ~100 bits** (`pcs-profile`, vultr3):

| N | scheme | proof | commit | open | verify |
|---|---|---|---|---|---|
| 2¹⁶ | BaseFold | 2.02 MB | 1.38 s | 1.25 s | 118 ms |
| 2¹⁶ | **Quasar** | 2.97 MB | **0.36 s** | **0.115 s** | 198 ms |
| 2²⁰ | BaseFold | 2.98 MB | 23.3 s | 20.2 s | 170 ms |
| 2²⁰ | **Quasar** | 4.43 MB | **5.74 s** | **1.85 s** | 280 ms |

So in this fork Quasar's win is **prover time against BaseFold — 3.9× commit and 10.8× open at
2²⁰** — while costing ~1.5× the proof size and ~1.6× the verify. That is the same shape as the
paper's claim (13.6–20.6× commit, 5.2–63.8× eval), scaled down because our BaseFold baseline is a
rate-1/4 random foldable code over ed25519 rather than their rate-1/2 setting.

**End-to-end, Fibonacci 2⁶ × 2 over ed25519** (stacked arity 7 — far below any crossover):

| backend | proof | prove | verify |
|---|---|---|---|
| Brakedown | 13.9 KB | 6.1 ms | 6.4 ms |
| Quasar | 125 KB | 12.6 ms | 44 ms |

Quasar's fixed code-switching cost (two BaseFold openings) dominates at toy sizes, exactly as the
paper's Fig. 2 shows its proof-size advantage over Brakedown only appearing from ≈2²¹. **The open
question E4 must answer: where is the crossover for the shapes this fork actually proves?** The
estimator can locate it without running the large sizes — that is what it is for.

### E5 — extend the cost model past the PCS

`quasar_cost_model.py` models the PCS alone. Add the multi-stark terms (zerocheck rounds, the layout
sumcheck, the per-table opening evaluations) so it predicts whole-proof bytes and verifier cost,
validated the same way against `pcs-profile`-style counters on the eddsa proof.

### E6 — optional: recursion readiness

The OpenVM recursion line needs every constant to be host/zkVM portable. Quasar's QA code and
BaseFold code are already ChaCha-seeded; the residual is the pre-existing ed25519 Poseidon2 constants
in the default hash path (`--features bb-hash` avoids them). `--features quasar,bb-hash` is then the
recursion operating point: measure guest instructions against the Brakedown and BaseFold baselines,
where Quasar's polylog verifier is the reason to care.

---

## 4. Build & test on vultr3

`ssh vultr3` — vast.ai container, 24 vCPU, 15 GB RAM, 30 GB free, RTX 3060 (8 GB), Rust 1.98 stable
at `~/.cargo/bin` (not on the non-interactive `PATH`). `/root/alp` already holds an older snapshot.

- Sync: `rsync -az --delete --exclude target --exclude .git alp2/ vultr3:/root/alp2/`.
- Build off the sshfs mount, always: `CARGO_TARGET_DIR=/root/alp2-target`.
- vast.ai ssh is flaky: `scp` files rather than heredocs; run long builds detached
  (`setsid nohup … >log 2>&1 </dev/null &`) and poll the log; strip the banner with
  `grep -vaE "Welcome to vast|Have fun|AI agents"`.
- **RAM is the binding constraint**: ed25519 elements are 32 B, so `N = 2²⁵` is ~1 GB for `M` plus
  ~1 GB for `C` at `c = 2`, before the Merkle tree. Plan sweeps to 2²²–2²⁴ and treat 2²⁵ as the
  ceiling; the paper's 2²⁹ rows need the 512 GB box and are out of reach here.
- Threads: benchmark at 24 threads (`--features parallel`), and record the 1-thread number too so
  the estimator has a thread-free baseline to match.

---

## 5. Python estimator + validation

`quasar_cost_model.py`, in the same style as `basefold_cost_model.py` / `brakedown_cost_model.py`
(CONFIG dict, optional JSON param file, `--json` machine-readable output).

**Inputs:** `num_vars`, `log_rows` (n₁), `c`, `field_bits`, `sec_bits`, `dist_failure_bits`,
BaseFold `inv_rate`/`reps`, hash preset (`babybear16` / `ed25519_native`).

**Modelled quantities** (each one traceable to a line of the Rust, as the other two models are):

| phase | model |
|---|---|
| security | δ_QA (M1 formula), `t = ⌈λ/−log₂(1−δ/3)⌉`, BaseFold reps, union bound (paper Eq. 2) |
| commit | `N` muls, `2cN log N₂` adds, `4N−1` leaf/compress hashes (Table 1) → perms via the hash preset |
| eval prover | `2N` (row combinations) + `2c·N₂` muls + `4cN₂ log N₂` adds (two encodes) + BaseFold commits over `2c` length-`N₂` vectors + `2·n₂` sumcheck rounds + selector sumcheck over `cN₂` |
| proof | `t·N₁` column elements + `t·log(cN₂)` digests + `4·(2n₂)` sumcheck coeffs + `2c` commitments + one batched BaseFold proof (reuse the `basefold_cost_model` term) |
| verifier | `t·N₁` muls/adds (column folding) + `t·log(cN₂)` compressions + `t·N₁·perms_per_elem` leaf hashes + `O(n₂)` for `H̃`/`eq` + batched BaseFold verify |

**Measured primitives.** Rather than fitting the arithmetic terms, the per-primitive costs are
measured once by `cargo run -p pcs-profile --release --bin ed25519_primitives` and used as
constants: a Fermat inversion is 77 mul + 252 sq; one degree-2 sumcheck-round interpolation is
170 mul + **2** inversions; `basefold::interp_deg2` inverts `2` afresh every round (77 mul);
sampling one ed25519 challenge costs 9 muls of Horner recombination; and the split leaf hash costs
**zero** field multiplications (the limb split is pure bit-shifting).

**Validation loop** (this is deliverable #4):
1. `pcs-profile/src/quasar.rs` emits `PcsMetrics` — proof bytes, hash perms by role
   (leaf/compress/sponge), ed25519 `mul/sq/inv/add` (`--features count-ops`), and wall-clock
   (`--no-default-features`). The counters are exact, not sampled.
2. `cargo run -p pcs-profile --release --bin sweep -- quasar 12,16,20,22` on vultr3.
3. Diff table `estimate vs measured` per column, per `num_vars`.
4. **Acceptance:** counted quantities (`t`, field elements, hash codes, digests) must match
   **exactly** — they are combinatorial. Serialized proof bytes and verifier perms/ops within
   **≤ 3 %** (the other two models sit at 0.2–1 % after a single serde-framing constant).
   Wall-clock is *not* modelled; it goes in the status log as measured only.
5. Any gap > 3 % is a bug in the model *or* the Rust — chase it before recording the row. Record the
   final calibration constants inline in the script's docstring, as the other two models do.

Cross-check against the paper where shapes allow: with `field_bits = 127, c = 2, N₁ = 2⁸, λ_dist =
103` the estimator must reproduce Table 4's `δ_QA` and `t` (already true in M0), and its proof-size
curve should track Fig. 2's shape (≈ 8 MB at `k = 29`) even though we can't run that size here.

---

## 6. Risks

| risk | mitigation |
|---|---|
| `t ≈ 420` columns × `N₁` rows dominates the proof | keep `N₁` small (2⁵–2⁶) and let the shape chooser optimize; `c = 4` cuts `t` to ~256 |
| our `basefold/` crate can't batch (M6) | M6 is scoped as a first-class milestone, not an afterthought; multi-column MMCS keeps it mechanical |
| degree-3 sumcheck over ed25519 is mul-heavy | carry `H̃` as a weight instead of a materialized poly; measure both |
| QA distance bound is *conjecture-free but external* ([LLX⁺26b] Thm 4.4 / Cor 3.21) | we port their formula and pin it with the Table 4 fixture; the paper's own numbers are the only reference |
| memory ceiling on vultr3 (15 GB) | cap sweeps at 2²⁴; the estimator extrapolates beyond |
| artifact is anonymized (license `LICENSE_MIT` in the mirror) | we *reimplement against Plonky3 traits* and cite the artifact; no code is copied verbatim into the repo |

---

## 7. Decisions taken

- **Field: ed25519 only** (2026-09-06). No `mersenne-127` crate. Consequence: the estimator is
  cross-checked against the paper's *security* numbers (Table 4) but our wall-clock is not directly
  comparable to Fig. 2 — comparisons here are Quasar vs our own Brakedown/BaseFold/QAPCS, which is
  the comparison that matters for this fork anyway.
- **Crate layout: new `qa-code/` member** (2026-09-06), not a file inside `brakedown/`.
- Inverse rate `c` stays a parameter; benchmark both 2 and 4 (c=4 cuts `t` from ~425 to ~256 at the
  cost of a 2× longer codeword).

---

## Status log

**2026-09-06 — batched end-to-end comparison, 1 → 256 signatures.** `eddsa-example --bin measure`
takes a signature count now: `EdDsaAir` stacks `2^k` signatures into one table (its preprocessed
selectors already restart per instance), and `src/trace.rs` builds the stacked main trace from
`curve::valid_instances`. Every stacked instance is fully constrained; only instance 0 is bound to
the public values, which is the same boundary the single-signature path already has.

All four backends, `num_bits = 253`, BabyBear hash, **serial** (no `parallel` feature), vultr3 —
raw data in `bench_data/eddsa_batch_bb_hash.json`:

| sigs | trace | | Brakedown | Ligero | BaseFold | **Quasar** |
|---|---|---|---|---|---|---|
| **1** | 256×45 | proof | 3.78 MB | **0.71 MB** | 2.87 MB | 3.96 MB |
| | | prove | 0.24 s | 0.35 s | 0.81 s | **0.25 s** |
| **16** | 4096×45 | proof | 6.32 MB | **2.64 MB** | 4.51 MB | 6.45 MB |
| | | prove | 2.7 s | 14.3 s | 13.1 s | **3.7 s** |
| **64** | 16384×45 | proof | 9.34 MB | 5.21 MB | **5.47 MB** | 8.00 MB |
| | | prove | **10.3 s** | 104.3 s | 53.1 s | 15.2 s |
| **256** | 65536×45 | proof | 15.14 MB | 10.34 MB | **6.51 MB** | 9.70 MB |
| | | setup | 7.5 s | 83.9 s | 375.8 s | **16.4 s** |
| | | prove | **40.3 s** | 795.6 s | 215.6 s | 53.6 s |
| | | verify | 1146 ms | 1771 ms | **374 ms** | 598 ms |
| | | per sig | 59 KB / 158 ms | 40 KB / 3108 ms | **25 KB** / 842 ms | 38 KB / **209 ms** |

Reading it:

- **Against BaseFold — the scheme Quasar shares a verifier profile with — Quasar is 4.0× faster to
  prove and 22.9× faster to set up at 256 signatures**, for 1.49× the proof. That is the paper's
  claim reproduced on a real AIR.
- **Against Ligero, 14.8× faster to prove at 256**, and the gap widens with every doubling:
  Ligero's superlinear Vandermonde encode takes its cost from 346 ms/sig to 3108 ms/sig, while
  Quasar's stays flat (252 → 209 ms/sig).
- **Against Brakedown, Quasar's proof crosses below between 16 and 64 signatures** (ratio 1.02× →
  0.86× → 0.64×), and its verifier is 1.9× faster at 256. Brakedown still proves 1.3× faster.
- Batching amortises what you would expect: Quasar's *prover* is already flat per signature, so
  what batching buys is proof size (3.96 MB → 38 KB/sig) and verification (2.3 ms/sig at 256).

**Correction to the earlier crossover estimate.** The Fibonacci `crossover` bin put the
Quasar-vs-Brakedown proof-size crossover near 2²³. That measurement was unfair: the fib example's
Brakedown backend used a hard-coded 32-query demo count while Quasar ran at ~100 bits. With both at
distance-derived query counts on the real AIR, the crossover is **between 16 and 64 signatures —
around 2¹⁹–2²⁰ evaluations**, much earlier than that estimate and earlier than the paper's ≈2²¹.
The fib config has been fixed to derive its query count the same way.

**2026-09-06 — Phase 2 done (E1–E6): Quasar proves the real EdDSA STARK.**

- **E1.** `eddsa-example/src/config.rs` restructured: the field/hash/transcript stack moved to
  `hash.rs`, and the commitment scheme to `backend/{brakedown,ligero,basefold,quasar}.rs` — one
  module per PCS exporting `CodePcs` + `build_pcs` + `NAME`. `config.rs` is now backend-free
  (the `ProverData<Self>` projection normalises to the concrete type, so
  `committed_table` needs no per-backend arm). All four backends prove and verify
  `EdDsaAir`, preprocessed trace and all: `tests/eddsa_roundtrip.rs` (renamed, backend-agnostic).
- **E2.** `QuasarParams::choose` picks the row/column split minimising predicted proof bytes over
  the shapes whose codeword can hold the target query count; `choose_best_effort` answers below
  that threshold with a clamped count. Pinned to the Python model on 12 reference shapes
  (`chooser_matches_cost_model`). This mattered: the hand-picked fallback it replaced produced a
  763 KB proof at the fib demo size where the chooser gets **15 KB**.
- **E3.** `QuasarParams::achieved_security_bits` reports the soundness a shape *delivers*, not the
  target; both configs log it at build time. A demo-sized instance clamped to a 2⁷ codeword is far
  under target, and now says so.
- **E4.** `eddsa-example --bin measure` and `ed25519-fib-example --bin crossover`. Numbers below.
- **E5.** `quasar_cost_model.py --stark LOG_HEIGHT:WIDTH[,WIDTH]` sums one PCS instance per
  committed trace. The multi-stark framework's own contribution is under 1% of a PCS-dominated
  proof, so it is deliberately **not** modelled: measured error is +0.7…+0.9% at fib arity 11–18,
  −2.3% on the two-instance EdDSA proof, and −7.8% only where the whole proof is 15 KB.
- **E6.** `--features quasar,bb-hash` — the recursion operating point — works. Both public codes are
  already ChaCha-seeded; `bb-hash` removes the last non-portable constants (the SmallRng-seeded
  ed25519 Poseidon2 round constants).

**End-to-end EdDSA proof** (`num_bits=253`, trace 256×45 + 8 preprocessed, vultr3):

| backend | ed25519 hash |  |  | BabyBear hash (`bb-hash`) |  |  |
|---|---|---|---|---|---|---|
| | proof | prove | verify | proof | prove | verify |
| Brakedown | 3.34 MB | 612 ms | 1317 ms | 3.78 MB | 243 ms | 235 ms |
| Ligero | **0.70 MB** | 503 ms | 164 ms | **0.71 MB** | 348 ms | **59 ms** |
| BaseFold | 2.51 MB | 2245 ms | 956 ms | 2.87 MB | 809 ms | 170 ms |
| **Quasar** | 3.58 MB | **558 ms** | 1237 ms | 3.96 MB | **253 ms** | 258 ms |

Read honestly: **Quasar beats BaseFold — the scheme it shares a verifier profile with — by 4.0× on
prove and 9.6× on setup**, and matches Brakedown's prover. It does *not* beat Ligero or Brakedown on
proof size at this trace, and it was never going to: 2¹⁴ is far below the size where a
polylogarithmic proof overtakes a square-root one.

**Where the crossover is** (`crossover` bin, Fibonacci over ed25519, Quasar vs Brakedown proof size):

| stacked arity | 11 | 15 | 18 | → |
|---|---|---|---|---|
| quasar / brakedown | 32.0× | 14.7× | 5.6× | **≈1× near 2²³** |

The ratio falls by ≈2.6× per four arity steps (Brakedown grows as `√N`, Quasar as `log²N`), putting
the crossover near **2²³** in this configuration. The paper sees ≈2²¹ — consistent, and later here
because our field elements serialize to 33 bytes rather than 16 and our BaseFold baseline is rate
1/4 with 288 repetitions rather than rate 1/2 with 249.

**So: use Quasar where the trace is large and the prover is the bottleneck; Ligero still wins the
small-trace proof-size contest in this fork.** That is a real result, not a disappointment — it is
the first time this fork can make that comparison at all.

**2026-09-06 — E0 done: Quasar proves a real STARK end to end.** `ed25519-fib-example` now has a
Quasar backend (`config_for_quasar`) and two integration tests: an honest Fibonacci proof
round-trips through `multi-stark` + `QuasarMlPcs` over ed25519, and a violated transition is
rejected by the zerocheck. The spike found exactly one class of missing code — `Clone` on
`QuasarProverData` and `QuasarMlProverData`, which the multi-stark prover requires — and nothing
else; the transcript and the prescribed-point opening worked on the first run. Also added
`tests/pcs_comparison.rs`, which prints proof size and prove/verify time per backend on one trace.
Phase 2 (§3b) plans the rest: the eddsa preprocessed-trace config, parameter selection from the
cost model, honest security accounting, and the crossover measurement.

**2026-09-06 — cost model done and validated.** `quasar_cost_model.py` now carries both halves;
`python3 quasar_cost_model.py --check bench_data/quasar_measured.json` diffs every prediction
against the profiler. Over **16 configurations** (num_vars 12–20 × log_rows 4–8 × index 2/4),
**0 mismatches**:

| quantity | worst error |
|---|---|
| `t`, `hash_leaf`, `hash_compress`, `hash_sponge`, `inv`, `sq`, `mul` | **exact** |
| `add` | 0.30 % |
| `proof_bytes` | 0.99 % |

Seven of the nine quantities are exact because they are derived from the code rather than fitted —
including the sponge count, which comes from replaying the duplex state machine, and the verifier's
multiplications, which use per-primitive costs measured by the new `ed25519_primitives` bin. Only
two numbers are fitted: a residual in the verifier's additions, and the mean postcard byte cost of a
digest (BabyBear limbs are varint-encoded, so digest size is genuinely data-dependent). Half of the
configurations were held out and checked afterwards — the accuracy is the same out of sample.

Two things the estimator says that the measurements could not reach:

- **Shape.** Sweeping `N₁` at 2²⁵ and 2²⁹ gives a shallow optimum at `N₁ = 2⁶–2⁷` — the paper's own
  benchmark choice — and `c = 4` beats `c = 2` on proof size at large `N` (≈8.6 MB vs 9.1 MB at
  2²⁹) at the cost of a 2× longer codeword to encode.
- **Extrapolation.** Predicted proof size at 2²⁹ is 8.6 MB against the paper's measured 8.05 MB —
  the same ballpark, which is as much as can be claimed across a different field (33-byte ed25519
  elements vs 16-byte Mersenne127) and a different hash.

**2026-09-06 — M1–M7 done. Quasar verifies end-to-end over ed25519.**

| milestone | landed as | tests |
|---|---|---|
| M1 | `security/src/qa.rs` (`QaRegime`, distance bound, query count, the full Eq. 2 union bound) | 6 new; **paper Table 4 reproduced exactly — δ_QA, t, λ_dist, λ_samp, λ_fld and λ_total, all ten rows** |
| M2 | new `qa-code/` crate (`wht`, `wht_batched`, `QaCode: LinearCode + LinearCodeExt`, `generator_matrix`) | 13 (encoder = generator matrix = dyadic convolution; systematic; linear; batched = per-column; `encode_ext`) |
| M2b | QAPCS for free: `TensorPcs<QaCode>` over ed25519 | 5 (roundtrip at c=2 and c=4, tampered column, wrong value) |
| M3–M5, M7 | new `quasar/` crate: `params.rs`, `util.rs` (H̃/eq helpers), `sumchecks.rs` (`Deg2Prover`/`Deg3Prover` on `p3-sumcheck`'s generic-degree engine), `pcs.rs` | 10 unit + 13 integration |
| M6 | `basefold`: `commit_many` / `commit_many_silent` / `open_batch` / `verify_batch` / `verify_batch_opening` (one width-`2k` pair-matrix tree ⇒ one Merkle path per query for the whole batch); `quasar`: the claim-reduction sumcheck | 4 new basefold tests; all 8 pre-existing ones still pass |
| M7 | `QuasarMlPcs` (`MultilinearPcs` + `PrescribedPointPcs`, drop-in beside `BasefoldMlPcs`), `pcs-profile/src/quasar.rs`, `quasar/benches/quasar_pcs.rs` | 103 tests green across every touched crate and consumer (incl. `eddsa-example`) |

Two deviations from the reference implementation, both simplifications that keep
the protocol faithful:

- The proximity and evaluation branches are batched by **one** challenge `β`
  *before* the encoding sumchecks, and the blocks by `α`, so the whole QA
  encoding check is a single degree-3 sumcheck plus a single degree-2 one —
  rather than the reference's per-instance layout. One shared `γ` as in the
  paper (Alg. 3 step 9), not one per block.
- Every evaluation claim (`4c+1` online + `c−1` preprocessed) is collapsed by
  one degree-2 reduction sumcheck onto a single point `ζ`, then discharged by
  exactly **two** batched BaseFold openings (online blocks, preprocessed
  coefficients).

**M2 measured — QA encoding beats both existing codes over ed25519** (vultr3, 24
threads, 8 rows encoded together):

| message length | QA rate ½ | Brakedown | RS rate ½ |
|---|---|---|---|
| 2¹⁰ | **1.42 ms** | 3.49 ms (2.5×) | 92.6 ms (65×) |
| 2¹² | **3.97 ms** | 12.27 ms (3.1×) | 1.476 s (372×) |

The 2.5–3.1× over Brakedown sits inside the paper's claimed 2.2–6.0× commitment
advantage; RS is superlinear, as expected.

**M7 measured — the whole PCS** (`pcs-profile`, ed25519 data + BabyBear-Poseidon2
hash, `N₁ = 2⁶`, query counts from the distance bound, ~103 bits):

| N | rate | columns `t` | proof | commit | open | verify |
|---|---|---|---|---|---|---|
| 2¹⁶ | ½ | 413 | 2.97 MB | 355 ms | 115 ms | **198 ms** |
| 2¹⁸ | ½ | 417 | 3.65 MB | 1.43 s | 456 ms | **238 ms** |
| 2²⁰ | ½ | 421 | 4.43 MB | 5.78 s | 1.85 s | **280 ms** |
| 2²⁰ | ¼ | 255 | 4.10 MB | 11.65 s | 2.91 s | **257 ms** |

**Verify grows 198 → 280 ms while `N` grows 16×** — the polylogarithmic verifier
behaving as advertised. Proof size is dominated by the two BaseFold openings'
Merkle paths (~1 MB each at 288 reps) plus `t·N₁·32 B ≈ 860 KB` of opened
columns; commit is hashing-bound, as it is for the Brakedown path in this fork.

Next: the cost model (§5) and its validation against these counters.

---

**2026-09-06 — M0 done.**

**2026-09-06 — M0 done.** Paper read (construction §3, efficiency App. A, parameters App. C).
Authors' artifact located and fully mirrored to `../quasar-ref/` (124 files). QA distance/query
formulas ported to Python and **validated exactly against Table 4, all 10 rows** (k=20…29:
δ_QA 0.444951→0.409635, t 445→487). Field survey done: 31-bit fields inadmissible for QA, ed25519
best (δ ≈ 0.46 at N₂=2¹⁶, t = 425 @ c=2; δ ≈ 0.73, t = 256 @ c=4). vultr3 checked out as the build
box (24 vCPU / 15 GB / RTX 3060 / Rust 1.98).
