# Quasar PCS — end-to-end benchmark report

Batched ed25519 signature verification (1 → 256 signatures in one proof) over the native ed25519

> **Individual benchmark runs now live in [`bench_results/`](bench_results/)**, one record per run
> with full provenance and an explicit comparability statement; the clustered change log is
> [`OPTIMIZATIONS.md`](OPTIMIZATIONS.md). This document is the narrative arc.

base field, proved four ways: **Brakedown**, **Ligero**, **BaseFold**, and **Quasar**
([eprint 2026/1839](https://eprint.iacr.org/2026/1839)).

Measured 2026-09-06, re-measured after optimisation and parallelisation 2026-09-07. Raw data:
[`bench_data/eddsa_batch_parallel.json`](bench_data/eddsa_batch_parallel.json) (**current — 24
threads**), [`eddsa_batch_optimised.json`](bench_data/eddsa_batch_optimised.json) (serial, optimised)
and [`eddsa_batch_bb_hash.json`](bench_data/eddsa_batch_bb_hash.json) (serial, pre-optimisation).
Implementation plan and milestone history: [`QUASAR_PCS_PLAN.md`](QUASAR_PCS_PLAN.md).

---

## 1. Hardware

| | |
|---|---|
| host | `vultr3` — vast.ai container (unprivileged Docker, not a VM) |
| CPU | Intel Xeon E5-2673 v3 @ 2.40 GHz, **24 vCPU** |
| RAM | 15 GiB (13 GiB available) |
| GPU | NVIDIA RTX 3060, 8 GiB — **unused**, no CUDA path is built yet |
| disk | 30 GiB free (overlay) |
| OS | Ubuntu 24.04.4 LTS |
| toolchain | rustc 1.98.0 stable, `--release`, `CARGO_TARGET_DIR=/root/alp2-target` |

> **§5 and §8 are single-threaded**; **§9 is the parallel re-measurement** and is the one to quote.
> §5/§8 are kept because the serial numbers isolate algorithmic cost from scheduling.

---

## 2. What is being proved

The monolithic `EdDsaAir` from `p3-ed25519-air`: **non-cofactored ed25519 signature verification**,
`R + [k]A == [s]B`, with subgroup checks on `A` and `R`, over the **native ed25519 base field**
`p = 2²⁵⁵ − 19` (non-two-adic, so no FFT-based scheme applies). The multilinear
SuperSpartan-flavoured `multi-stark` prover drives it.

**Batching.** `EdDsaAir.log_instances = k` stacks `2^k` signatures vertically into one table. Its
preprocessed trace restarts the row selectors (`is_instance_first` / `is_instance_last`) at every
instance boundary, and the per-instance constraints are gated on those rather than on the global
first row, so **every stacked signature is fully constrained**. The stacked main trace is built by
[`eddsa-example/src/trace.rs`](eddsa-example/src/trace.rs) from `curve::valid_instances`, which
gives each instance *distinct* secrets — the trace is not `n` copies of one witness.

| signatures | trace rows | main columns | preprocessed columns | stacked arity (main) |
|---|---|---|---|---|
| 1 | 256 | 45 | 8 | 2¹⁴ |
| 4 | 1 024 | 45 | 8 | 2¹⁶ |
| 16 | 4 096 | 45 | 8 | 2¹⁸ |
| 64 | 16 384 | 45 | 8 | 2²⁰ |
| 256 | 65 536 | 45 | 8 | 2²² |

Two commitments per proof — the main trace and the preprocessed trace stack different column counts,
so each gets its own PCS instance sized for its own arity.

### Scope boundaries (stated, not hidden)

- **Only instance 0 is bound to the public values.** The AIR reads them `when_first_row`, so a
  batched proof asserts "I know 2^k valid signatures, the first of which is this one" rather than
  exposing 2^k public keys. Every other instance is still fully constrained — the test
  `eddsa_batched_roundtrips_and_rejects_a_corrupt_instance` corrupts the *second* instance and
  checks the proof is rejected, under all four backends.
- **`k = H(R‖A‖M)` is not proved.** The challenge scalar is a public input, not a hashed one; the
  SHA-512 half of RFC 8032 lives outside this AIR. Pre-existing, unrelated to the PCS comparison.
- **Poseidon2 round constants are demo-grade** on the ed25519 hash path (seeded from a fixed RNG).
  The BabyBear path used here uses OpenVM's fixed Horizen RC16 constants instead.

---

## 3. Hash function — the answer to "Poseidon2 over BabyBear?"

**Yes.** These runs used `--features bb-hash`, which is a *cross-field* stack: the committed data,
the codes and every sumcheck stay over **ed25519**, and only the hashing moves to **BabyBear**.

| role | construction |
|---|---|
| permutation | **Poseidon2 over BabyBear**, width 16, with OpenVM's fixed Horizen RC16 constants (`default_babybear_poseidon2_16`) — bit-identical to the OpenVM poseidon2 precompile |
| Merkle leaf hash | `SplitFieldHasher<Ed25519, BabyBear, PaddingFreeSponge<·, 16, 8, 8>>` — each ed25519 element is split into **9 BabyBear limbs** by radix-2³⁰ bit shifts (injective, no BigUint division), then absorbed at rate 8 |
| node compression | `TruncatedPermutation<·, 2, 8, 16>` — one permutation per Merkle node |
| digest | 8 BabyBear ≈ 248 bits ⇒ ~124-bit collision resistance |
| transcript | `DuplexChallengerCF<Ed25519, BabyBear, ·, 16, 8>` — sponges over BabyBear, presents an ed25519 interface (9-limb split on observe, Horner recombination on sample) |

Why: a native 255-bit Poseidon2 permutation costs hundreds of emulated 64-bit-limb multiplications
inside a 32-bit zkVM, while a BabyBear permutation is native `u32` arithmetic — ~100× cheaper per
permutation, and a candidate precompile. Opened leaves stay ed25519 values, so no reconstruction is
needed on the arithmetic side.

The alternative (default, no `bb-hash`) is **Poseidon2 over ed25519**, width 3, rate 2, 1-element
digest (~127-bit collisions). It is kept for host-side comparison; it makes proofs ~10% smaller
(33-byte digests instead of 39) and roughly 2–5× slower.

---

## 4. Code parameters

### Shared

| | |
|---|---|
| base field | ed25519, `p = 2²⁵⁵ − 19` (255 bits) |
| challenge field | the same field (degree-1 self-extension; `\|F\| ≈ 2²⁵⁵` is ample) |
| scalar bit length | `num_bits = 253` (the real ed25519 size) |
| layout | `PrefixProver` stacked witness, padding floor 1 |
| serialization | postcard; an ed25519 element is 33 bytes, a BabyBear digest ≈ 39 |

### Per backend

| backend | code | rate | queries | shape rule |
|---|---|---|---|---|
| **Brakedown** | recursive Spielman (Haböck `small_proof` preset, δ = 0.070) | ~1/1.42 | `t = ⌈100 / −log₂(1−δ/3)⌉ = 2936`, clamped | `params::optimal_log_rows` |
| **Ligero** | Reed–Solomon | 1/2 | `t = ⌈100 / log₂ rate⌉ = 100`, clamped (list-decoding-capacity **conjecture**) | square split |
| **BaseFold** | random foldable | 1/4 | 288 (list-decoding radius δ ≈ 0.215, 100-bit target) | n/a — commits the stacked poly directly |
| **Quasar** | QA (quasi-abelian), index `c = 2` | 1/2 | `t` from the proven QA distance at **λ_dist = λ_samp = 103** ⇒ 411–427 columns, clamped | `QuasarParams::choose` — minimises predicted proof bytes |
| ↳ code switching | random foldable | 1/4 | 288 | message length `N₂` |

Quasar's column count is derived, not chosen: `δ_QA` comes from the [LLX⁺26b] bound (Thm 4.4 /
Cor 3.21) as implemented in [`security/src/qa.rs`](security/src/qa.rs), which reproduces the paper's
Table 4 exactly. Over ed25519 at these sizes `δ_QA ≈ 0.46`, giving `t ≈ 420` at rate 1/2. `N₁` (the
matrix row count, and hence the height of an opened column) is picked by `QuasarParams::choose`,
pinned to the Python cost model on 12 reference shapes.

Both public codes are **ChaCha-seeded**, so a prover on x86 and a verifier on riscv32 build
identical parameters — a requirement for the OpenVM recursion path.

### Reproducing

```bash
# The whole sweep, including the JSON parser that regenerates bench_data/:
CARGO_TARGET_DIR=/local/path bench_data/run_batch_sweep.sh

# Or one backend directly (PROFILE=1 also prints the prover's span tree):
cargo run -p eddsa-example --release --bin measure \
  --features quasar,bb-hash,parallel -- 253 1,4,16,64,256
```

---

## 5. Results

All values: `num_bits = 253`, BabyBear hash, **serial**, vultr3.

### Proof size (MB)

| signatures | Brakedown | Ligero | BaseFold | **Quasar** |
|---|---|---|---|---|
| 1 | 3.78 | **0.71** | 2.87 | 3.96 |
| 4 | 4.71 | **1.35** | 3.65 | 5.09 |
| 16 | 6.32 | **2.64** | 4.51 | 6.45 |
| 64 | 9.34 | 5.21 | **5.47** | 8.00 |
| 256 | 15.14 | 10.34 | **6.51** | 9.70 |

### Prover wall-clock (s)

| signatures | Brakedown | Ligero | BaseFold | **Quasar** |
|---|---|---|---|---|
| 1 | **0.24** | 0.35 | 0.81 | 0.25 |
| 4 | **0.76** | 2.10 | 3.24 | 0.94 |
| 16 | **2.68** | 14.27 | 13.08 | 3.74 |
| 64 | **10.27** | 104.27 | 53.06 | 15.21 |
| 256 | **40.33** | 795.59 | 215.60 | 53.60 |

### At 256 signatures — the full picture

| | Brakedown | Ligero | BaseFold | **Quasar** |
|---|---|---|---|---|
| setup | 7.5 s | 83.9 s | 375.8 s | **16.4 s** |
| prove | **40.3 s** | 795.6 s | 215.6 s | 53.6 s |
| verify | 1146 ms | 1771 ms | **374 ms** | 598 ms |
| proof | 15.14 MB | 10.34 MB | **6.51 MB** | 9.70 MB |
| **per signature** | 59 KB / 158 ms | 40 KB / 3108 ms | **25 KB** / 842 ms | 38 KB / **209 ms** |

### How each scheme scales

| | prove per signature (1 → 256 sigs) | reading |
|---|---|---|
| Brakedown | 236 → 158 ms | flat, linear-time code |
| **Quasar** | 252 → **209 ms** | flat, linear-time code |
| BaseFold | 810 → 842 ms | flat but ~4× more expensive |
| Ligero | 346 → **3108 ms** | **superlinear** — Reed–Solomon's Vandermonde encode |

Quasar's proof crosses **below** Brakedown's between 16 and 64 signatures (ratio 1.02× → 0.86× →
0.64×), i.e. around 2¹⁹–2²⁰ evaluations.

> **Correction to an earlier estimate.** A Fibonacci-based sweep put that crossover near 2²³. It was
> unfair: the fib example's Brakedown backend used a hard-coded 32-query demo count while Quasar ran
> at ~100 bits. With both at distance-derived counts on the real AIR, the crossover is much earlier.
> The fib config has been fixed.

---

## 6. Analysis — which backend for which bottleneck

**Quasar wins the prover contest outright among the polylog-verifier schemes.** Against BaseFold —
the scheme it shares a verifier profile with — it is **4.0× faster to prove and 22.9× faster to set
up** at 256 signatures, which is the paper's central claim reproduced on a real AIR. Against Ligero
it is **14.8× faster**, and the gap widens with every doubling because Ligero's encode is
superlinear while Quasar's is linear.

**But BaseFold still has the smaller proof and the cheaper verifier** at 256 signatures (6.51 MB /
374 ms vs 9.70 MB / 598 ms), and **Brakedown still proves 1.3× faster**. So the choice depends on
which cost actually binds:

| if the bottleneck is… | pick | why |
|---|---|---|
| prover throughput / cost per signature | **Quasar** | 4× BaseFold, 15× Ligero; flat 209 ms/sig; and it overtakes Brakedown on both proof size and verify beyond ~64 signatures |
| proving latency for one small batch | Brakedown | fastest prover below the crossover, simplest stack |
| smallest proof at small batches | Ligero | 0.71 MB at one signature — but note its query count rests on a **conjecture**, and its prover collapses when batched |
| recursion guest cost today | BaseFold | smallest proof (the guest's input) and cheapest verifier at 256 sigs |

**For a strategy built on batching many signatures and recursing:** Quasar is the right *prover*,
and it is close enough on the verifier side that the remaining gap is an engineering matter rather
than a structural one — §7 shows the proof can plausibly be cut by ~70% with changes that touch
parameters rather than the protocol, which would put it below BaseFold on both axes. Today, if the
guest-side verifier is the hard constraint and prover time is free, BaseFold is still ahead by a
nose. **If prover time is not free — and at 216 s vs 54 s per 256-signature batch it rarely is —
Quasar is the better choice.**

One caveat before acting on this: nothing here has been measured *inside* a zkVM guest. The
recursion metric that matters is guest instructions/cells, not host milliseconds, and the mapping is
not linear (hash permutations dominate, and a poseidon2 precompile would change the ranking).
`--features quasar,bb-hash` is exactly the configuration to measure that with, and it works today.

---

## 7. Optimisation roadmap

### Where the proof size actually goes

For the 256-signature main instance (2²² stacked), from the validated cost model:

```
proof 5.27 MB   =   1.16 MB  field elements  (22%)
                +   4.12 MB  Merkle digests  (78%)   ← of which 93% are BaseFold query paths
```

**93% of the digests — over 70% of the whole proof — are Merkle authentication paths from the two
BaseFold openings.** Every lever below attacks that term.

| lever | effect on the 2²² instance | cost | effort |
|---|---|---|---|
| **Merkle cap height 8** (`MerkleTreeMmcs::new(_, _, 8)`, currently 0) | **−51%** → 2.56 MB | commitment carries 256 digests instead of 1 | one constant |
| **Quasar's *inner* foldable-code rate, 1/4 → 1/16** (bigger radius ⇒ 288 → 161 queries) — the `BF_INV_RATE` constant in `backend/quasar.rs`, **not** the standalone BaseFold backend | **−30%** → 3.72 MB | 4× longer inner codeword, but on length-`N₂` vectors — see the note below | one constant |
| **24-bit grinding** (the sumcheck machinery already takes `pow_bits`) | **−24%** → 4.01 MB | ~16 M hashes once, per proof | wire `pow_bits` through |
| **QA index `c = 4`** (δ 0.46 → 0.73 ⇒ t 425 → 256) | −7% | 2× longer QA codeword, ~2× commit | one constant |
| **all four together** | **−72% → 1.47 MB** | | |

> **Two different BaseFolds — worth separating.** Quasar *contains* BaseFold: it commits its `2c`
> length-`N₂` blocks with a foldable code and discharges every claim through two batched BaseFold
> openings. That inner instance is what the second row above retunes. It is a different object from
> the **BaseFold backend** the results compare against, which commits the whole stacked polynomial.
>
> This matters because the rate lever is **asymmetric in Quasar's favour**. At the 2²² instance the
> chooser picks `N₁ = 2⁶`, so Quasar's inner code works on `N₂ = 2¹⁶` vectors while the standalone
> backend works on all 2²²:
>
> | foldable-code work (elements encoded) | rate 1/4 | rate 1/16 |
> |---|---|---|
> | Quasar's inner BaseFold (4 blocks × `N₂`) | 1.0 M | **4.2 M** |
> | standalone BaseFold backend (full 2²²) | 16.8 M | 67.1 M |
>
> **Quasar at rate 1/16 still encodes 4× less than the BaseFold backend does at rate 1/4** — and the
> backend at 1/16 would need a 2.1 GB codeword before the Merkle tree. Code switching is precisely
> what makes an aggressive foldable-code rate affordable: it moves the per-element cost onto a
> vector 64× shorter.
>
> The other three levers (Merkle cap, grinding, and — for the tensor schemes — the code rate) apply
> to **every** backend and would shrink all of them, so they do not by themselves change the ranking.
> Only the inner-rate lever is structurally Quasar's.

Extrapolated to the whole 256-signature proof: **9.70 MB → ≈2.7 MB**, which would put Quasar below
BaseFold's 6.51 MB *and* below Ligero, while keeping its prover advantage. These are parameter
changes, not protocol changes — all four are already expressible in the existing code.

A fifth, structural one: the preprocessed QA coefficient vectors are a **separate** BaseFold
commitment needing a second full opening, and that opening costs as much as the online one despite
committing 8× less data. Giving it its own (much higher) rate — it is preprocessed, so encode time
is paid once, ever — would nearly halve the switching cost on its own.

### Prover speed

| lever | expected | effort |
|---|---|---|
| **`--features parallel`** | large — 23 of 24 cores are currently idle, and the Merkle build (the dominant commit cost) is already rayon-ready | one flag; **do this first** |
| **Parallelise the Quasar sumcheck round provers** | 2–4× on the sumcheck phase. `Deg2Prover`/`Deg3Prover::interpolate` allocates a `Vec` per factor per round and runs serially; fold in place and compute `h(0)`, `h(2)`, `h(3)` in one pass | half a day |
| **Hoist `basefold::interp_deg2`'s `TWO.try_inverse()`** | small but free — it runs a **full Fermat inversion (77 muls + 252 squarings) every round** just to compute ½ | ten minutes |
| **Carry `H̃(γ,·)` as a tensor weight** instead of materialising a length-`N₂` table in the outer sumcheck | one fewer full pass + `O(N₂)` memory | half a day |
| **Cache the preprocessed commitment** in the verifying key instead of recomputing it in `QuasarPcs::new` | most of the 16.4 s setup | small |
| **GPU WHT encoder** (the paper reports 33–148×; the artifact ships the CUDA kernel, and vultr3 has an RTX 3060) | limited until leaf hashing also moves to the GPU — commit is Merkle-bound in this configuration, not encode-bound | days |

### Suggested order

1. `--features parallel` and re-run this whole sweep — it changes every absolute number here.
2. Merkle cap 8 (all backends, for fairness) + Quasar's inner foldable rate 1/16 — two constants,
   ~62% off Quasar's proof, and the second one is a lever only Quasar can afford.
3. Measure the recursion metric that actually decides the strategy: guest instructions for the
   Quasar verifier under `--features quasar,bb-hash`, against the BaseFold and Brakedown baselines.
4. Then decide whether grinding and `c = 4` are worth their prover cost.

---

## 8. Optimisation round 1 — pruning, caps, grinding

Three of the levers from §7 are now implemented. Two paid off; one did not, which is worth as much
as the two that did.

### What was done

| lever | where | applies to |
|---|---|---|
| **Merkle path pruning** | `basefold/src/merkle_pcs.rs` (one pruned multi-opening per fold round) and `brakedown/src/pcs.rs` (one for the queried columns), both via the generic `Mmcs::{open_multi_batch, verify_multi_batch}` | **all four** backends |
| **Grinding** (16 bits) | `BasefoldPcs::new_with_grinding` + `QuasarPcs::with_grinding`; ground after the final value is absorbed and *before* the query indices are drawn, so a prover who dislikes its query set must redo the work | BaseFold and Quasar |
| **Merkle caps** | `hash::MERKLE_CAP_HEIGHT` | tried at 6, **reverted to 0** |

### Why caps were reverted

A cap looks like free money against *unpruned* paths. Against pruned ones it is a **net loss**:
pruning already removes the top levels — with hundreds of queries the upper ~log₂(queries) levels
are covered by some other query and their siblings are recomputable — so the cap buys nothing off
the paths while adding `2^h` digests to **every** commitment, and Quasar has many of those (column
tree, block tree, coefficient tree, and one fold-oracle root per round per opening).

Measured at cap 6: the EdDSA proof *grew* from 2.61 → 2.69 MB at 16 signatures and 3.33 → 3.43 MB at
64. The constant is kept at 0 with the finding recorded next to it.

### Results — before → after

| backend | sigs | proof | verify | prove |
|---|---|---|---|---|
| Brakedown | 1 | 3.78 → **1.01 MB** (−73%) | 216 → **94 ms** (−56%) | 0.24 → 0.24 s |
| | 256 | 15.14 → **12.28 MB** (−19%) | 1146 → **1058 ms** (−8%) | 40.3 → 41.1 s |
| Ligero | 1 | 0.71 → **0.65 MB** (−8%) | 59 → **45 ms** (−24%) | 0.35 → 0.35 s |
| | 256 | 10.34 → **10.28 MB** (−1%) | 1771 → **1742 ms** (−2%) | 795.6 → 799.5 s |
| BaseFold | 1 | 2.87 → **0.78 MB** (−73%) | 171 → **52 ms** (−70%) | 0.81 → 1.52 s |
| | 256 | 6.51 → **2.50 MB** (−62%) | 374 → **160 ms** (−57%) | 215.6 → 216.7 s |
| **Quasar** | 1 | 3.96 → **1.32 MB** (−67%) | 258 → **80 ms** (−69%) | 0.25 → 1.38 s |
| | 256 | 9.70 → **3.82 MB** (−61%) | 598 → **271 ms** (−55%) | 53.6 → 56.2 s |

Pruning is worth **−61% proof and −55% verify** for Quasar at 256 signatures, and it helps exactly
in proportion to how path-dominated a proof was: Quasar and BaseFold (−61…−73%) are almost all
paths, Brakedown less so (−19% at scale, its proof being mostly opened column *data*), and Ligero
barely at all (−1%) since its 100 queries share almost no ancestors.

**Grinding costs a fixed ~1.1 s** (one grind per query phase, several phases per proof) for a ~14%
proof reduction. That is a bad trade at one signature (+446% prove time) and a good one at 256
(+5% prove for −14% proof) — it amortises, so it suits the batching strategy but should be turned
off for single-signature latency. `POW_BITS = 20` was measurably too slow under the 255-bit ed25519
hash and was reduced to 16.

### The state of play at 256 signatures

| | Brakedown | Ligero | BaseFold | **Quasar** |
|---|---|---|---|---|
| setup | 7.7 s | 84.1 s | 376.2 s | **16.5 s** |
| prove | **41.1 s** | 799.5 s | 216.7 s | 56.2 s |
| verify | 1058 ms | 1742 ms | **160 ms** | 271 ms |
| proof | 12.28 MB | 10.28 MB | **2.50 MB** | 3.82 MB |
| per sig | 48 KB / 161 ms | 40 KB / 3123 ms | **10 KB** / 847 ms | 15 KB / **219 ms** |

**Quasar now beats Brakedown on proof size from 4 signatures onward** (ratio 1.31 → 0.94 → 0.64 →
0.44 → **0.31**), where before optimisation the crossover was between 16 and 64. Against Brakedown at
256 signatures it is **3.2× smaller and 3.9× faster to verify**, for 1.37× the prover time. Against
BaseFold it remains **3.9× faster to prove and 22.7× faster to set up**, at 1.53× the proof and
1.7× the verify.

### Fairness note

Pruning was applied to **all four** backends. Grinding was wired into BaseFold and Quasar only —
they share the machinery — so Brakedown and Ligero are *under*-optimised here by roughly the ~14%
that grinding buys. Plumbing proof-of-work through `TensorPcs` would close that gap and is the
honest next step before quoting these ratios anywhere load-bearing.

### What is left from §7

- **`--features parallel`** — still the biggest single lever, still untouched. 23 of 24 cores idle.
- **Quasar's inner foldable-code rate 1/4 → 1/16** — the one lever that is structurally Quasar's
  (see the note in §7), now worth relatively more since pruning has removed the shared paths and
  what remains scales with the query count.
- QA index `c = 4`; the sumcheck-prover and `interp_deg2` micro-optimisations; GPU WHT.

---

## 9. Parallel execution (current numbers)

Everything above ran on one core. Making it use 24 took two kinds of work: flipping the rayon shim,
and rewriting the loops that were written serially — plus one algorithmic fix that turned out to
matter more than all the threading.

### The algorithmic fix, first

`RandomFoldableCode::new` inverted its twist vector **one element at a time**, and there are
`inv_rate · 2^log_len` of them. Over ed25519 an inversion is `a^(p−2)` — 252 squarings and 77
multiplications — so building the code for a 2¹⁹-variable BaseFold instance was ~2 M Fermat
inversions. That, not encoding, was BaseFold's entire setup cost.

Replacing it with `batch_multiplicative_inverse` (Montgomery's trick: **one** inversion plus ~3n
multiplications for the whole level) took BaseFold's setup at 16 signatures from **22.2 s to 0.81 s**
— 27×, and none of it from threads. Quasar's inner instance benefits too, on shorter vectors.

### What was parallelised

| crate | loop |
|---|---|
| `basefold` (had **no** rayon dependency at all) | the foldable-code encode recursion (`join` above level 8) and fold; `eq_table`; `bind_msb`; the per-round `h(0)/h(1)/h(2)` reduction, fused into one pass; the eta-combination; the pair-matrix interleave |
| `quasar` | the sumcheck round provers, rewritten as a **single parallel reduction per round** that computes every reported node at once (the old shape interpolated each factor at each node into a fresh `Vec`); `tensor_table`; `combine_rows`; `claim_weights`; the elementwise combines |
| features | `basefold` gained a `parallel` feature; `eddsa-example`/`ed25519-fib-example` now forward it to every backend explicitly — `p3-basefold`'s was **not** reachable from `p3-multi-stark`, so it would have stayed serial |

### Where the time actually went — profiling, not guessing

An earlier draft blamed the modest 5.8× on "grinding, 11%". Amdahl's law does not support that: if
1.1 s of a 56.2 s serial prover were the only serial part, 24 cores would give
`55.1/24 + 1.1 ≈ 3.4 s`, not 9.8. So the phases were instrumented
(`PROFILE=1 cargo run -p eddsa-example --bin measure …` prints the span tree) and the container's
CPU budget checked (`cpu.max` = 2304000/100000 = **23.04 cores**, so the cores are real).

Measured breakdown of a 4.71 s prove at 64 signatures:

| phase | time | share |
|---|---|---|
| **proof of work — six separate grinds** | **~2.6 s** | **~55%** |
| column-tree Merkle build (2 instances) | 0.54 s | 11% |
| multi-stark zerocheck sumcheck | 0.57 s | 12% |
| QA encode + block commits | 0.42 s | 9% |
| everything else (row combine, 4 sumchecks, openings) | ~0.6 s | 13% |

A Quasar proof runs **six** independent grinds — a column phase plus two BaseFold opening phases, for
each of the two committed traces — and `GrindingChallenger::grind`'s default implementation tests
candidates in a **serial loop**. At 55% serial, Amdahl caps the speedup at ~1.8× no matter how many
cores exist. That, not scheduling overhead, was the answer.

### Parallel proof of work

`basefold::grind_parallel` searches windows of candidates in parallel and returns the **smallest**
witness, not the first one found — so proofs stay bit-identical across thread counts and across the
two rayon backends. It leaves the transcript in exactly the state the serial `grind` would, which is
required rather than cosmetic: the verifier's `check_witness` performs the same absorption.

Quasar's prove at 256 signatures fell **9.77 → 6.96 s**; at 64 signatures **4.68 → 2.23 s** (2.1×).
BaseFold gained little (its two grinds were a small share of a much larger prover).

### Serial → parallel, at 256 signatures

| backend | setup | prove | speedup |
|---|---|---|---|
| Brakedown | 7.7 → **3.53 s** | 41.1 → **7.44 s** | 5.5× |
| Ligero | 84.1 → **10.42 s** | 799.5 → **77.29 s** | 10.3× |
| BaseFold | 376.2 → **10.86 s** (34.6×) | 216.7 → **22.43 s** | 9.7× |
| **Quasar** | 16.5 → **3.32 s** | 56.2 → **6.96 s** | **8.1×** |

### The table that matters

256 signatures, 24 threads, BabyBear hash, all optimisations on:

| | Brakedown | Ligero | BaseFold | **Quasar** |
|---|---|---|---|---|
| setup | 3.53 s | 10.42 s | 10.86 s | **3.32 s** |
| prove | 7.44 s | 77.29 s | 22.43 s | **6.96 s** |
| verify | 1076 ms | 900 ms | **164 ms** | 271 ms |
| proof | 12.28 MB | 10.28 MB | **2.50 MB** | 3.82 MB |
| **per signature** | 48 KB / 29 ms | 40 KB / 302 ms | **10 KB** / 88 ms | 15 KB / **27 ms** |

**Quasar now wins the prover outright** — faster than Brakedown (6.96 vs 7.44 s) while producing a
**3.2× smaller proof** and verifying **4.0× faster**, and **3.2× faster than BaseFold**, the only
other polylogarithmic-verifier scheme here. It is also the fastest to set up. On every axis except
raw proof size (where BaseFold's 2.50 MB still leads) it is the best option in the tree.

### What is left

Quasar's 8.1× on ~23 cores is 35% efficiency, and the remaining loss is now identifiable rather than
mysterious:

- **Merkle build and zerocheck** are the two largest phases and both already use rayon; their
  scaling is bounded by memory bandwidth on this 2.4 GHz Xeon, not by missing parallelism.
- **Sumcheck round tails** operate on 2–8 elements, where a task costs more than the work. A serial
  cutoff below some length would trim the overhead.
- **Six grinds is itself the question worth asking** — the security accounting needs each query phase
  gated, but 16 bits × 6 phases is a design choice, not a requirement. Fewer bits, or fewer phases,
  is a cheaper answer than making the search faster.

Then, still untouched from §7: Quasar's inner foldable-code rate 1/4 → 1/16 (the structurally-Quasar
proof-size lever), QA index `c = 4`, and the GPU WHT encoder.

---

## 10. Recursion: real OpenVM compression of a Quasar proof

Everything above measures the **native** proof. This section closes the loop the whole line of
work exists for: hand a Quasar proof to an OpenVM guest, have the guest verify it, and prove
*that* — collapsing the proof into a single constant-size recursive proof.

### Setup

| | |
|---|---|
| Box | vultr4 — AMD EPYC 7B12 (30.72-core cgroup quota), 503 GB RAM, **RTX 5090 32 GB** (sm_120, driver 580.159.03), CUDA 13.1 |
| Toolchain | `cargo-openvm v2.0.2 (59a69b8) [cuda]`, guest `riscv32im-risc0-zkvm-elf`, nightly-2026-01-18 |
| Guest | `alp-openvm/guests/quasar-verifier/` — rebuilds the whole `--features "quasar,bb-hash"` stack in-circuit |
| Host blob | `eddsa-example --bin gen_proof` → `(num_bits, log_instances, log_height, pis, vk, proof)` via postcard |
| Precompiles | rv32i/rv32m/io + **modular arithmetic** (p3-ed25519 routes field mul through `moduli_declare!`) |

Nothing but the proof, the verifying key and the public values crosses the boundary. Both codes
(QA and foldable), the Poseidon2 permutation and the `QuasarParams` shape are ChaCha- or
constant-seeded, so the guest reconstructs them and any divergence surfaces as a transcript
mismatch rather than a silent pass.

### Results

| | 1 signature | 4 signatures |
|---|---|---|
| Quasar proof in | 1,320,413 B | 1,776,130 B |
| `log_height` | 8 | 10 |
| guest verdict | `[1, …]` verified | `[1, …]` verified |
| guest instructions | 1,192,414,643 | 2,081,714,400 |
| segments | 84 | 144 |
| `prove stark` (GPU) | 202 s | 316 s |
| **aggregated proof out** | **630,732 B** | **630,732 B** |
| `verify stark` | ✅ | ✅ |

The output is the *same 630,732 bytes* in both cases — the aggregated proof's size is fixed by the
recursion circuit, not by what it verified. Compression is therefore 2.1× at one signature, 2.8× at
four, and about **6.1× at 256 signatures**, where the native Quasar proof is 3.82 MB (§9). Guest
cost over the same 1→4 range grows only 1.75×, well under the 4× in signatures.

### Where the guest spends its cycles

`guests/quasar-verifier/prof/` is the same guest with a second input word selecting how far to run,
so the counts subtract cleanly:

| phase | instructions | share |
|---|---|---|
| postcard deserialize of the 1.32 MB blob | 92,218,545 | 7.7 % |
| **rebuild both PCS instances** | **248,568,134** | **20.8 %** |
| `multi_stark::verify` proper | 851,628,440 | 71.4 % |

**A fifth of the guest is setup that does not belong there.** `QuasarPcs::verify_core` never reads
`self.code`, and the BaseFold verifier needs only `diag_at(level, j)` at the queried positions — yet
the guest materialises both codes in full. `RandomFoldableCode::new` samples `inv_rate·(2^log_cols −
1)` ed25519 twists and batch-inverts them; `QaCode::new_random` samples `(c−1)·2^log_cols`. Worse,
this term *grows with the batch*: at 256 signatures `log_cols` reaches ~15, so the eager path would
sample ~131k field elements before verifying anything.

Two fixes, both invisible to the host:

1. a verifier-only `QuasarPcs` constructor that takes the coefficient commitment directly instead of
   re-deriving it;
2. a lazy `RandomFoldableCode` verifier path (`new_verifier` / `diag_at`) deriving one twist on
   demand from the ChaCha stream.

### Two obstacles worth recording

- **One OpenVM version per guest graph.** `p3-ed25519` pulls `openvm` at **tag v2.0.1** under
  `[target.'cfg(target_os = "zkvm")'.dependencies]` for the modular precompile. A guest asking for
  v2.0.2 gets both compiled, each defining `#[panic_handler] fn panic_impl`, and the link fails with
  `E0152: found duplicate lang item 'panic_impl' … first defined in crate std` — which points at
  serde and std, not at the actual cause. Pin the guest to v2.0.1.
- **`VerifyingKey` was not serializable.** `multi-stark/src/keys.rs` now derives `Clone, Serialize,
  Deserialize` with `#[serde(bound = "")]`, so the guest receives the preprocessed commitment rather
  than regenerating and re-committing the preprocessed trace in-circuit.

### Not yet done

`cargo openvm prove evm` — the final wrap into a halo2 SNARK of a few hundred bytes — needs
`cargo openvm setup --evm` (large download plus a heavy keygen). The STARK aggregation above is the
step that establishes recursion works; the EVM wrap is a packaging layer on top of it.

---

## 11. Recursion round 2 — removing the guest's `O(N)` setup, and Quasar vs BaseFold in-circuit

§10 measured a Quasar proof verifying in OpenVM and found a fifth of the guest was setup rather
than verification. This section closes that, fixes an unfair baseline, and answers the question
that actually decides the strategy: **which PCS is cheapest inside the recursion circuit?**

### Correcting the earlier baseline

The 33 s / 163 M-instruction figure quoted in §10's discussion came from `gen_proof **16**` — a
16-bit toy scalar used as a plumbing smoke test. At the real 253-bit size the same Ligero guest
costs **621,108,478 instructions**, matching the 621 M recorded independently in
`alp-openvm/chips/poseidon2/integration/README.md`. Nothing regressed; the comparison was against
a smaller problem.

### The Poseidon2 precompile was not enabled

The guest's `openvm.toml` declared only `rv32i/rv32m/io/modular`. The chip is built and validated
on this box but was not wired in. Adding `[app_vm_config.poseidon2]` plus a
`Poseidon2Precompile` permutation (same Horizen RC16 constants as the host's
`default_babybear_poseidon2_16`, so transcripts still agree) is worth ~1.9× to both guests:

| guest | no chip | with chip |
|---|---|---|
| Ligero | 621,108,478 instr / 110 s | 328,557,786 / 53 s |
| Quasar | 1,192,414,643 / 189 s | 629,558,262 / 106 s |

### Two `O(N)` setup terms, and only one of them was worth removing

**(a) The coefficient recommit — removed, a clear win.** `QuasarPcs::new` ends with
`basefold.commit_many_silent(code.scaled_coeffs())`, a full encode-and-Merkle-commit of the public
coefficient vectors, of which the verifier reads only the digest. `QuasarPcs::new_verifier` now
takes that digest instead, `QuasarMlPcs::new_verifier` drops the prover key, and the QA code
becomes `Option` — a verifier instance never samples it. The host ships the two digests beside the
verifying key (they are public parameters fixed by the code seed). Test
`verifier_only_instance_matches_the_full_one` pins accept-and-reject equivalence.

**(b) The foldable-code twists — implemented, then measured as a loss for Quasar.**
`RandomFoldableCode::new` samples `inv_rate·(2^n − 1)` twists; the verifier reads `O(q·log n)` of
them. Sampling is now *addressable* (a ChaCha stream per element, seeded `seed + idx·φ`) and
`new_verifier` derives twists on demand via `diag_at`. `lazy_matches_full` pins the two codes equal.
For plain **BaseFold** this is the difference between an `O(N)` and an `O(q log N)` guest. For
**Quasar** it is a net *loss*: its inner code is over `N₂ = 2^log_cols`, only 2^9 at one signature
and 2^11 at sixteen, while the query count is a size-independent 241 — so `q·log N₂ ≈ 2169`
positions are touched to avoid materialising 2044, each paying a Fermat inversion. Measured at
S=1: **557 M eager → 822 M lazy**. The Quasar guest therefore keeps the eager code; the finding is
recorded in the guest source so it is not "fixed" again.

### Head-to-head in the recursion circuit

`num_bits=253`, Poseidon2 + modular precompiles, identical hash stack, grinding and blob format;
Quasar with the verifier-only constructor and eager inner code, BaseFold with the lazy code.

| | Quasar S=1 | BaseFold S=1 | Quasar S=16 | BaseFold S=16 |
|---|---|---|---|---|
| host proof in | 1,320,300 B | 781,776 B | 2,243,487 B | 1,485,539 B |
| guest instructions | 559,940,571 | **482,961,447** | 1,262,729,474 | **761,412,765** |
| segments | 48 | **41** | 107 | **64** |

**BaseFold is the cheaper scheme to verify in-circuit, and its lead widens with the batch** — 1.16×
at one signature, 1.66× at sixteen. This inverts the host ranking, where Quasar proves 3.2× faster
than BaseFold (§9). The reason is structural: Quasar's opening is *two* BaseFold openings plus a
column phase, and its 348-of-1024 column sample is set by the security target rather than by the
trace size, so at recursion-sized witnesses it pays a large constant that its polylog asymptotics
have not yet earned back.

### What this means for the strategy

Quasar and BaseFold answer different questions, and the honest split is:

- **Native proving of a large batch** — Quasar. Fastest prover in the tree, fastest setup, proof
  smaller than Brakedown's from 4 signatures on.
- **The recursion / aggregation layer** — BaseFold, by a widening margin.

They compose: nothing stops a Quasar proof of the witness from being aggregated under a BaseFold
recursion layer, since the guest rebuilds only public parameters. The levers that would change the
verdict are the ones that cut Quasar's *query counts* — a larger QA index `c`, an inner BaseFold
rate of 1/16, or more grinding bits — not further guest engineering.

---

## 12. Making the recursion guest cheap (2026-09-10)

§11 left the Quasar guest at 646.9 M instructions for 16 signatures and concluded BaseFold was the
better recursion scheme. **That conclusion is now stale**: the same guest is 87.1 M, a **7.4×**
cut, and none of the five changes touched the scheme. Every one was found by measuring where the
instructions went, and three of them contradicted a plausible-sounding prediction.

### The five

| # | change | where | effect at 16 sigs |
|---|---|---|---|
| 1 | `QuasarParams::choose` moved to the host and shipped in the blob | `backend/quasar.rs`, `gen_proof` | −95 M |
| 2 | Public parameters baked into the program image at build time | `guests/quasar-verifier/build.rs`, `RandomFoldableCode::new_verifier_with_tables` | −309 M → 2 M |
| 3 | Deserialiser's `BigUint` round-trips removed | `Ed25519::{from,to}_canonical_bytes` | −92 M |
| 4 | Struct given the modular precompile's alignment | `#[cfg_attr(target_os = "zkvm", repr(C, align(32)))]` | −5.4 M |
| 5 | Blob moved to a flat, borrowed-bytes format | new crate `flat-serde` | −17.4 M |

**(1) An `f64` bisection in a circuit with no FPU.** `QuasarParams::choose` bisects on `log`, `pow`
and `expm1` to pick four integers. On rv32 every one is a softfloat call. The shape is a public
parameter of a fixed trace geometry, so the host resolves it and ships it beside the coefficient
commitments.

**(2) Public parameters are public — bake them.** The foldable code's twists are a deterministic
function of a seed, yet a zkVM has no cross-run state, so every proof re-derived all ~65k of them.
Sampling one costs 6,809 instructions; reading one from a `static` costs **34**. `build.rs` emits
the table using the same `sample_nonzero_at` the verifier would have called, and
`baked_sampler_matches_the_derived_code` pins them equal.

*Bake the inverses too.* With only the twists baked the guest got **worse** (177 M → 368 M at one
signature): sampling became free, and the lazy path then paid a 14,879-instruction field inversion
per queried position. `TINV_WORDS` is what makes lazy strictly better than eager rather than a
trade — and it retires the eager-vs-lazy question that had been muddy for two rounds of
measurement.

**(3) The same disease, one layer down.** `impl Deserialize for Ed25519` allocated a `Vec<u8>`,
built two `BigUint`s, **rebuilt the modulus as a `BigUint` on every call**, and round-tripped
through bytes again — ~1,830 instructions to read 32 bytes. The replacement does four `u64` reads
and a four-word `< P` comparison. No wire-format change.

**(4) 34 instructions of chip, 144 of copying.** `Ed25519Fp` is `#[repr(C, align(32))]`; a
`[u64; 4]` is the same 32 bytes but 8-aligned, so every operation *copied* into the precompile's
type and back. Giving the struct the chip's alignment makes `zk::as_fp` a pointer cast. Microbench:
multiply 178 → 66, square 186 → **34**, the bare chip op.

**(5) Parsing, not the format.** postcard decodes a `Vec<Ed25519>` at 710 instructions per element
against 152 for a flat read — and the gap is serde's per-element machinery, not varints, so a
codec swap would not have captured it. `p3-flat-serde` keeps the derived structure but hands the
visitor **borrowed** bytes, so a field element is parsed in place. The blob grows ~5 %.

### Where it lands

| signatures | instructions | segments | `prove app` | app proof | per signature |
|---|---|---|---|---|---|
| 1 | 54,882,084 | 6 | — | — | — |
| 16 | 87,103,545 | 9 | 40 s | 6.9 MB | 2.50 s |
| 32 | 96,352,089 | 10 | — | — | — |
| 64 | 106,067,766 | 11 | 45 s | 8.5 MB | 0.70 s |
| 128 | 116,577,012 | 13 | 45 s | 9.5 MB | **0.35 s** |

About **+12 M instructions per doubling** — genuinely logarithmic, so per-signature cost keeps
falling. Native proving for comparison at 128 signatures: 2.66 s, 20.8 ms/sig, 3.24 MB.

### Three predictions that were wrong

Worth recording, because the failures were more instructive than the successes.

- **A Merkle-path chip.** Argued for at "~80 % of the guest", then 26 %, then measured at
  **8.5 %** — the first figure came from `verify_batch`, which the guest never calls (it uses the
  pruned `verify_batch_pruned`), and the second from a per-path cost at a depth most rounds do not
  reach. Not built, and it should not be.
- **The coefficient opening costing 3.4× the online one.** An ablation artifact: the two are the
  same function called twice, and dead-coding the second let the compiler specialise the first. A
  counter showed identical work (24 calls / 5,784 query-paths each). **Ablation is only valid when
  a cut removes code that exists nowhere else.**
- **Op-count × microbench-cost.** Predicted ~36 M for change (4); it delivered 5.4 M. A
  `black_box`ed loop forces spills real code elides. Useful for ranking, never for promising.

### Still open

- **The preprocessed opening is 79 % the size of the main one**, for a trace that is fixed and
  public, and every ablation delta is a sum over both PCS instances — so it is roughly half the
  guest. Its row selectors are periodic and could be evaluated in closed form instead of opened;
  its doubling-ladder columns are arbitrary curve coordinates and cannot. Uninvestigated, and
  plausibly worth more than anything above.
- Deserialisation is still ~35 % of what remains, against a `read_vec` floor of 7.3 M.
- The inner code rate 1/4 → 1/16 would cut the inner query count 241 → 136, but lengthens the
  inner codeword 4× and so slows the *native* prover — declined, since that is Quasar's advantage.
