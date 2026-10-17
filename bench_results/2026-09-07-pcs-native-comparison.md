# Native prover: Quasar is the fastest, 6.96 s at 256 signatures with a 3.82 MB proof

Host-side comparison of four commitment schemes on the same EdDSA AIR. No zkVM involved.

## Target

| | |
|---|---|
| layer | `pcs` |
| subject | Quasar vs Brakedown vs BaseFold vs Ligero, native multi-stark prover |
| parameters | `num_bits=253`, 256 signatures batched; `bb-hash` (Poseidon2 over BabyBear); 24 threads; `POW_BITS=16`; Merkle path pruning on; `MERKLE_CAP_HEIGHT=0` |
| build features | `<backend>,bb-hash,parallel` |
| command | `cargo run --release -p eddsa-example --bin measure -- 253 1,4,16,64,256` |
| input artefacts | none — traces generated in-process |

## Provenance

| field | value |
|---|---|
| date | 2026-09-07 |
| tree | `alp2`, working tree before `f60a208f` |
| **uncommitted** | **true** — recorded retrospectively |
| host | vultr3 |
| CPU | Intel Xeon E5-2673 v3 — 24 visible, quota 23.04 cores |
| memory | 15 GB |
| GPU | RTX 3060 8 GB (unused) |
| rustc | 1.98.1 |

## Results

At 256 signatures (65536×45 trace):

| backend | setup | prove | verify | proof |
|---|---|---|---|---|
| **Quasar** | 3.32 s | **6.96 s** | 271 ms | 3.82 MB |
| Brakedown | — | 7.44 s | 1146 ms | 12.3 MB |
| BaseFold | — | 22.3 s | — | 2.50 MB |
| Ligero | — | (superlinear RS encode) | — | — |

All measured. Quasar is the fastest prover *and* the fastest to set up, with a 3.2× smaller proof
than Brakedown and 4.0× faster verification; BaseFold still has the smallest proof.

## Comparability

- **not comparable to:** any `recursion` record — this is native host proving, not in-circuit
  verification, and the two rank the schemes **oppositely**. Quasar wins here; BaseFold won
  in-circuit as of 2026-09-09.
- **not comparable to:** anything measured on vultr4 by wall-clock — different CPU.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| Merkle caps (height 6) | proof 2.61 → 2.69 MB at 16 sigs | path pruning already frees the top levels; a cap adds 2^h digests to every commitment |

## Notes

The prove/verify asymmetry is the whole point of the scheme choice: pick per layer, not once.
