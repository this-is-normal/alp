# Pricing the lazy code path: inversion is 77.9% of it, sampling 22.1% — and the hybrid fallback beats batch inversion

<!-- layer: recursion -->
<!-- target: Quasar recursion guest, foldable-code derivation cost -->

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | the cost of *deriving* rather than reading `RandomFoldableCode` elements in-guest |
| **parameters** | `num_bits=253`, **256 signatures**, `log_cols=16`, `BF_INV_RATE=4`, `POW_BITS=16`, `MERKLE_CAP_HEIGHT=0`, `MAX_LOG_COLS=16` |
| build features | guest `probe-derive-inv` / `probe-derive-both` (**cost attribution only**) · host blob `quasar,bb-hash,parallel` |
| command | `cargo openvm build [--features probe-…]` then `cargo openvm run --exe <vmexe> --input qz_in_256.json --mode segment` |
| input artefacts | `/root/eddsa_quasar_256.bin` -> `qz_in_256.json`, unchanged across all three variants |

**Why S=256 and not the cliff.** At 512/1024 the fallback is all-or-nothing: sampling and inversion
move together and cannot be separated. At 256 every index fits the baked table, so each operation
can be forced back on independently. Per-position cost does not depend on level, so the ratio
transfers to the cliff — and the per-level totals below confirm it does.

## Provenance

| field | value |
|---|---|
| date | 2026-09-11 |
| tree | `alp2` @ `70cb077d` (`quazar2`), **19 uncommitted files** |
| host | vultr4 — AMD EPYC 7B12 (quota 30.72 cores), 503 GB, RTX 5090, CUDA 13.1 |
| rustc | 1.98.1 |
| cargo-openvm | v2.0.2 (59a69b8) [cuda] |

## Results

All three variants reveal `[1, 0, …]` — the probes compute identical values by a costlier route.

| variant | code lookups | instructions | segments | `.vmexe` bytes |
|---|---|---|---|---|
| `base` | table twists + table inverses | **89,069,686** | 9 | 93,036,034 |
| `inv` | table twists + **derived** inverses | 272,944,784 | 26 | 49,920,250 |
| `both` | **derived** twists + **derived** inverses | 325,216,364 | 30 | 7,025,360 |

**`base` reproduces the independently measured 89,069,686 exactly**, so the probe scaffolding
(a new `build_code` function and three `cfg` branches) is cost-neutral and the deltas are
attributable to derivation alone.

| operation | instructions | share | per code level (16) |
|---|---|---|---|
| **inversion** (`inv − base`) | 183,875,098 | **77.9 %** | 11.49 M |
| **sampling** (`both − inv`) | 52,271,580 | **22.1 %** | 3.27 M |
| total derivation (`both − base`) | 236,146,678 | 100 % | **14.76 M** |

### Cross-check against the cliff

`2026-09-11-recursion-quasar-direct-preprocessed` derived **~14.4 M per level** from the 512/1024
excess over trend. This probe gives **14.76 M per level** at S=256 by an unrelated route. Agreement
within 2.5 % across two independent measurements.

### What the fixes are now worth

*(modelled from the measured per-level cost; the underlying per-level figure is measured)*

| fix | 1024 sigs | vs 367.4 M measured |
|---|---|---|
| batch inversion alone | ~167 M | 2.2x |
| **per-index hybrid table fallback alone** | **~140 M** | **2.6x** |
| both together | ~116 M | 3.2x |

**The hybrid wins alone and is the simpler build.** It removes *both* costs at 16 of 18 levels;
batch inversion shrinks only the 77.9 % share, everywhere. They compose, but if one gets built it
should be the hybrid.

### Table cost in image bytes, measured directly

| image | bytes | table data it carries |
|---|---|---|
| guest with no tables (`both`) | 7,025,360 | — |
| + twists (`inv`) | 49,920,250 | 8.4 MB |
| + twists and inverses (`base`) | 93,036,034 | 16.8 MB |

**~43 MB of image per 8.4 MB of field data — 5.1x inflation.** Scaling to `MAX_LOG_COLS=18`
(4x the data, 67 MB) gives ~343 MB + 7 MB base ≈ **350 MB**, confirming from measurement the figure
previously extrapolated from a single `MAX_LOG_COLS` 14->16 data point.

## Comparability

- **comparable to:** `2026-09-11-recursion-quasar-direct-preprocessed` — same box, commit,
  parameters and input blob; `base` here is that record's 256-signature row, reproduced exactly.
- **not comparable to:** any run at a different `log_cols` — the per-level totals scale with the
  number of levels, which is what the per-level column normalises away.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| deriving the per-position cost from microbench op costs (6,809 sample / 14,879 inverse) | inconsistent | implies more inversions than samples, which is impossible since every inverted position is first sampled. Confirms the documented rule that op-count x microbench-cost is an **upper bound**, not an estimate |
| measuring the split at 512/1024 | not attempted | the fallback there is all-or-nothing; the two operations cannot be separated |

## Notes

**Two earlier estimates of this split were wrong, in opposite directions.** First "~95 % is
inversion", then a correction to "~69 %". Measured: **77.9 %**. The correction was closer and its
qualitative claim (inversion dominant, sampling not negligible) was right, but neither number was
derived from measurement and neither should have been stated as if it were.

**The probe features are cost-attribution only** and are declared as such in the guest manifest,
mirroring alp2's `ablate-*` flags. Unlike those, they do not make verification return early — they
compute the same values by a more expensive route, which is why `[1, …]` is a valid check on them.

**Codegen-artifact risk was considered and is low.** `cfg`-selected branches are the shape that
produced a phantom 3.4x asymmetry in an earlier ablation. Here the variants differ only in which
`fn` pointers are stored in `RandomFoldableCode` and called indirectly from `p3-basefold`, and the
`base` control reproduces its independent measurement to the instruction.
