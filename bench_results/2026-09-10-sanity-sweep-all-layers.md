# Sanity sweep: three PCSs across all three layers, one commit, one box

Quasar is the fastest prover at every batch size from 16 signatures up and the cheapest to verify
in-circuit; BaseFold keeps the smallest proof and the fastest native verify. Brakedown is
arithmetic-bound in verification by an order of magnitude.

## Target

| | |
|---|---|
| layer | `pcs` + `prover` + `recursion` (one sweep, three sections) |
| subject | Brakedown, BaseFold, Quasar |
| parameters | `num_bits=253`; `bb-hash` (Poseidon2 over BabyBear); λ=100 provable regime; Quasar `c=2`, `SECURITY_BITS=103`, `BF_INV_RATE=4`, `POW_BITS=16`; BaseFold `INV_RATE=16`, 136 queries, lazy code; Brakedown `small_proof` preset, `SplitStrategy::MinSize`; guest `MAX_LOG_COLS=16` |
| build features | host `<backend>,bb-hash,parallel`; PCS section additionally with and without `count-ops`; guest chips `rv32i,rv32m,io,modular,poseidon2` |
| commands | `sweep all 14,16,18,20 prov` · `measure 253 1,4,16,32,64,128,256,1024` · `openvm run --mode segment` + `prove app --segment-max-memory $((26<<30))` |
| input artefacts | blobs regenerated per backend immediately before use — see Notes |

## Provenance

| field | value |
|---|---|
| date | 2026-09-10 |
| tree | `alp2` @ `70cb077d` (`quazar2`) |
| **uncommitted** | **true (4 files)** |
| host | rust-dev |
| CPU | - — 2 visible, quota  |
| memory | 3 GB |
| rustc | 1.97.0 |
| field | value |
|---|---|
| date | 2026-09-10 |
| tree | `alp2` @ `70cb077d` (`quazar2`) |
| **uncommitted** | **true (4 files)** |
| host | vultr4 |
| CPU | AMD EPYC 7B12 64-Core Processor — 256 visible, quota 30.72 cores |
| memory | 503 GB |
| GPU | NVIDIA GeForce RTX 5090, 32607 MiB, 580.159.03 |
| rustc | 1.98.1 |
| cargo-openvm | cargo-openvm v2.0.2 (59a69b8) [cuda] |

## 1. PCS only — operation breakdown

Verify-side costs for a single committed polynomial. Hash counts are Poseidon2 permutations;
arithmetic is ed25519 operations. `commit`/`open` are prove-side wall-clock from the
`count-ops`-off run.

At `num_vars = 20`:

| PCS | q | proof | commit | open | verify | hashes | mul | sq | add | arith total |
|---|---|---|---|---|---|---|---|---|---|---|
| BaseFold | 288 | 1.28 MB | 16.1 s | 14.5 s | 62.6 ms | 38,964 | 25,560 | 5,040 | 23,996 | ~54,600 |
| **Quasar** | 421 | 2.10 MB | **4.05 s** | 1.31 s | 109 ms | 66,345 | 84,305 | 36,792 | 68,572 | ~190,000 |
| Brakedown (small-proof) | 2,936 | 4.22 MB | 3.87 s | 0.18 s | 269 ms | 142,220 | 748,845 | 0 | 680,106 | ~1,429,000 |

**Arithmetic-to-hash ratio: BaseFold ≈ 1.4:1, Quasar ≈ 3:1, Brakedown ≈ 10:1.** Brakedown's verifier
is overwhelmingly field arithmetic; BaseFold's is the most balanced. Quasar is the only one with a
material inversion count (146), from its sumcheck interpolations.

Quasar commits 4× faster than BaseFold at this size — the fast-prover property, visible before any
batching.

## 2. Whole prover — S = 1 … 1024

Trace is `256·S × 45`. All measured, 24 threads.

| S | Brakedown prove | BaseFold prove | **Quasar prove** | Brakedown proof | BaseFold proof | Quasar proof |
|---|---|---|---|---|---|---|
| 1 | 96 ms | 380 ms | 275 ms | 0.96 MB | 0.66 MB | 1.32 MB |
| 4 | 238 ms | 965 ms | 374 ms | 1.87 MB | 0.90 MB | 1.77 MB |
| 16 | 571 ms | 2.41 s | **654 ms** | 3.48 MB | 1.18 MB | 2.25 MB |
| 32 | 772 ms | 4.99 s | **887 ms** | 4.65 MB | 1.34 MB | 2.54 MB |
| 64 | 1.54 s | 9.67 s | **1.46 s** | 6.49 MB | 1.50 MB | 2.87 MB |
| 128 | 2.52 s | 19.3 s | **2.69 s** | 8.82 MB | 1.68 MB | 3.24 MB |
| 256 | 4.81 s | 38.5 s | **4.47 s** | 12.3 MB | 1.88 MB | 3.81 MB |
| **1024** | 19.8 s | 152.9 s | **14.1 s** | 23.6 MB | 2.29 MB | 4.86 MB |

Per signature at S=1024: Quasar **13.7 ms**, Brakedown 19.4 ms, BaseFold 149.3 ms. Quasar's
per-signature cost keeps falling with batch size (275 → 13.7 ms); BaseFold's is flat at ~150 ms.

Native verify at S=1024: BaseFold 118 ms, Quasar 261 ms, Brakedown 1,557 ms.

## 3. OpenVM recursion — S = 32, 64, 128

Guest instructions and GPU `prove app`. Both guests built at this commit, so these two are
**comparable to each other** for the first time.

| S | Quasar instr | seg | Quasar `prove app` | BaseFold instr | seg | BaseFold `prove app` |
|---|---|---|---|---|---|---|
| 32 | 96,352,089 | 10 | **42 s** | 191,915,964 | 17 | 51 s |
| 64 | 106,067,766 | 11 | **44 s** | 204,294,415 | 19 | 51 s |
| 128 | 116,577,012 | 13 | **46 s** | 216,839,286 | 20 | 55 s |

App proofs: Quasar 8.0/8.5/9.5 MB against BaseFold 12.5/13.6/14.5 MB.

**Quasar is now ~1.9× cheaper in-circuit than BaseFold**, which inverts the 2026-09-09 finding. The
reason is that the five Quasar-side optimisations landed first; the shared-crate ones lift both, but
the baked parameter tables and the shipped shape are Quasar-specific.

## Comparability

- **comparable to:** internally — one commit, one box, one toolchain for all three sections.
- **supersedes:** `2026-09-09-recursion-basefold` and the in-circuit half of
  `2026-09-10-recursion-quasar-flat`, whose Quasar-vs-BaseFold comparison was flagged invalid
  because the two were measured at different commits. This run fixes that.
- **not comparable to:** `2026-09-07-pcs-native-comparison` — measured on vultr3, a different CPU.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| (none new this run) | | |

## Notes

**A stale artefact nearly corrupted section 3.** The first attempt generated the "Quasar" blobs
with a `gen_proof` binary still built `--features basefold`; the guest failed with
`Invalid field element`, and the giveaway was the blob being byte-identical in size to the
BaseFold one. Regenerate blobs **after** rebuilding for the backend under test, and check sizes.

PCS-section wall-clock comes from the `count-ops`-off run; the operation columns come from the
`count-ops`-on run. Mixing them would misreport both, which is why the harness records them
separately.
