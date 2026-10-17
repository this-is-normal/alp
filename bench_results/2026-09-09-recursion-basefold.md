# BaseFold recursion guest: 168.8 M instructions at 16 signatures, 0.227 s/signature at 256


> **Superseded for the Quasar-vs-BaseFold comparison** by [`2026-09-10-sanity-sweep-all-layers`](2026-09-10-sanity-sweep-all-layers.md), which measures both guests at one commit on one box. The in-circuit ranking there is the opposite of the one implied here.

**Superseded for comparison purposes** — see the Comparability section before using these numbers
against Quasar's.

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | BaseFold PCS over a random foldable code, verified in an OpenVM guest |
| parameters | `num_bits=253`; signatures 1/16/256; `INV_RATE=16`, `RADIUS=0.350` ⇒ 136 queries; `POW_BITS=16`; lazy foldable code (`new_verifier`) |
| build features | host `basefold,bb-hash,parallel`; guest chips `rv32i,rv32m,io,modular,poseidon2` |
| command | `gen_proof 253 <n>` then `cargo openvm prove app --input … --segment-max-memory $((26<<30))` |
| input artefacts | blobs from `gen_proof` at the same tree state |

## Provenance

| field | value |
|---|---|
| date | 2026-09-09 |
| tree | `alp2` @ working tree before `f60a208f` |
| **uncommitted** | **true (~25 files)** — recorded retrospectively; the tree was not committed at measurement time |
| host | vultr4 |
| CPU | AMD EPYC 7B12 — 256 visible, quota 30.72 cores |
| memory | 503 GB |
| GPU | NVIDIA GeForce RTX 5090, 32607 MiB, 580.159.03 |
| CUDA | 13.1 |
| rustc | 1.98.1 |
| cargo-openvm | v2.0.2 (59a69b8) [cuda] |

## Results

| signatures | instructions | segments | `prove app` | app proof | per signature |
|---|---|---|---|---|---|
| 1 | 122,800,769 | 11 | 36 s | 8.7 MB | 36 s |
| 16 | 168,826,803 | 15 | 50 s | 12.3 MB | 3.1 s |
| 256 | 219,601,831 | 19 | 58 s | 15.7 MB | **0.227 s** |

All measured; every run reveals `[1, …]`.

Phase split at 16 signatures, by ablation: deserialize 37.9 M (22 %), build config 0.3 M,
sumcheck replay 4.0 M, grinding + query sampling + Merkle 18.5 M (11 %), fold-chain walk
107.9 M (64 %). That 64 % was the lazy code's per-position field inversion.

## Comparability

- **not comparable to:** `2026-09-10-recursion-quasar-flat` — these predate three fixes that live
  in **shared** crates (`p3-flat-serde`, the `Ed25519` alignment, and the `BigUint`-free serde) and
  would lift BaseFold too. Any Quasar-vs-BaseFold in-circuit ratio drawn from these two records is
  wrong. Re-measure both at one commit.
- **not comparable to:** the vultr2-era BaseFold recursion figures (76/111/131 segments) — different
  box, and its ~1.09 s/segment against this box's ~2 s/segment reflects single-thread CPU speed,
  not the scheme.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| rate 1/4 instead of 1/16 | 241 queries instead of 136 | 1.8× more verifier work; the 4× longer codeword costs the prover little |
| eager foldable code | `O(N)` construction | loses badly at scale; the lazy path was the fix here (the opposite of Quasar's verdict) |

## Notes

BaseFold looked like the better recursion scheme when this was measured. That conclusion did not
survive the Quasar work of 2026-09-10 and should not be quoted from here.
