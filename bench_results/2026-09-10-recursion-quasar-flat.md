# Quasar recursion guest: 87.1 M instructions at 16 signatures, 0.35 s/signature at 128


> **Superseded for the Quasar-vs-BaseFold comparison** by [`2026-09-10-sanity-sweep-all-layers`](2026-09-10-sanity-sweep-all-layers.md), which measures both guests at one commit on one box. The in-circuit ranking there is the opposite of the one implied here.

Verifying a Quasar-committed Plonky3 EdDSA proof inside an OpenVM guest, after the day's five
optimisations. Down from 646.9 M instructions at the same batch size — **7.4×**.

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | Quasar PCS (QA code + code switching into BaseFold), verified in an OpenVM guest |
| parameters | `num_bits=253`; signatures 1/16/32/64/128; QA index `c=2`; `SECURITY_BITS=103`; `BF_INV_RATE=4`, `BF_RADIUS=0.215`; `POW_BITS=16`; `MERKLE_CAP_HEIGHT=0`; guest `MAX_LOG_COLS=16` |
| chooser output | main PCS `N₂ = 2^13/2^14/2^15/2^16` at 16/32/64/128 sigs; preprocessed two levels smaller |
| build features | host `quasar,bb-hash,parallel`; guest chips `rv32i,rv32m,io,modular,poseidon2` |
| command | `gen_proof 253 <n>` then `cargo openvm prove app --input … --segment-max-memory $((26<<30))` |
| input artefacts | blobs regenerated at this commit — the flat format and the baked table both change them |

## Provenance

| field | value |
|---|---|
| date | 2026-09-10 |
| tree | `alp2` @ `70cb077d` (`quazar2`) |
| **uncommitted** | **true (1 files)** |
| host | vultr4 |
| CPU | AMD EPYC 7B12 64-Core Processor — 256 visible, quota 30.72 cores |
| memory | 503 GB |
| GPU | NVIDIA GeForce RTX 5090, 32607 MiB, 580.159.03 |
| rustc | 1.98.1 |
| cargo-openvm | cargo-openvm v2.0.2 (59a69b8) [cuda] |

## Results

| signatures | instructions | segments | `prove app` | app proof | per signature |
|---|---|---|---|---|---|
| 1 | 54,882,084 | 6 | — | — | — |
| 16 | 87,103,545 | 9 | 40 s | 6.9 MB | 2.50 s |
| 32 | 96,352,089 | 10 | — | — | — |
| 64 | 106,067,766 | 11 | 45 s | 8.5 MB | 0.70 s |
| 128 | 116,577,012 | 13 | 45 s | 9.5 MB | **0.35 s** |

All measured. Every run reveals `[1, …]`; `verify app` passes. Growth is ~+12 M instructions per
doubling — logarithmic, so per-signature cost keeps falling.

Phase split at 16 signatures, by ablation: deserialize 36.5 M (33 %), build config 0.06 M,
column + selector + encoding sumchecks ~18 M (16 %), the two inner BaseFold openings 55.4 M (50 %).

### What produced the 7.4×

| change | effect at 16 sigs |
|---|---|
| `QuasarParams::choose` moved to the host, shape shipped in the blob | −95 M |
| public parameters baked into the program image (twists **and** inverses) | 309 M → 2 M |
| `Ed25519` serde without `BigUint` | −92 M |
| `#[repr(C, align(32))]` on `Ed25519` so `zk::as_fp` is a pointer cast | −5.4 M |
| blob moved to `p3-flat-serde` (borrowed bytes) | −17.4 M |

## Comparability

- **comparable to:** each other — one guest binary, one box, only the batch size differs.
- **not comparable to:** `2026-09-09-recursion-basefold` — those numbers predate the last three
  changes above, which live in shared crates and would lift BaseFold too. Re-measure before
  quoting a Quasar-vs-BaseFold in-circuit ratio.
- **not comparable to:** any pre-2026-09-09 recursion figure — a 33 s result quoted earlier was
  `gen_proof **16**`, a 16-bit toy scalar, not 253.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| Quasar inner rate 1/4 → 1/16 | 1.26 B → 2.13 B instructions | 4× longer inner codeword ⇒ 4× the eager twist sampling. Also slows the *native* prover, which is Quasar's advantage |
| lazy foldable code with per-position inversion | 557 M → 822 M at 1 sig | a 14,879-instruction field inversion per queried position beat the sampling it saved |
| baking twists but not inverses | 177 M → 368 M at 1 sig | sampling became free, the per-position inversion remained |
| `MAX_LOG_COLS=14` | 64 sigs: 128 M → 349 M | a shape above the table falls back to deriving, silently |
| `MAX_LOG_COLS` sized larger than needed | 16 sigs: 28 s → 44 s | a 92 MB `.vmexe` costs ~16 s per prove over a 27.7 MB one, on identical work |

## Notes

`MAX_LOG_COLS` is a tuning knob, not a constant: size it to the deployment's maximum batch — 14
for ≤32 signatures, 16 for 128.

Native proving for context at 128 signatures: setup 1.43 s, prove 2.66 s (20.8 ms/sig), verify
210 ms, proof 3.24 MB.
