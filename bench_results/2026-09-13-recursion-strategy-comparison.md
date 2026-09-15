# Three strategies for the baked twist table, compared: **C+D wins**; Merkle buys a 13x smaller image and loses on everything else

<!-- layer: recursion -->
<!-- target: Quasar recursion guest, twist-table strategies C / C+D / M -->

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | how the guest obtains the foldable code's twists and inverses |
| **parameters** | `num_bits=253`, S=32..1024, `BF_INV_RATE=4`, `POW_BITS=16`, `MERKLE_CAP_HEIGHT=0`, `MAX_LOG_COLS=16`, `MAX_ROOT_LOG_COLS=20` |
| commands | `cargo openvm run --exe … --mode segment` / `--mode meter`; `cargo openvm prove app --exe … --segment-max-memory $((26<<30))` |
| input artefacts | `c2_*` (C), `cd2_*`/`hint_quasar_*` (C+D), `merkle_*` (M) — all regenerated 2026-09-13 in the current blob format |

## The three strategies

| | twists | inverses | image |
|---|---|---|---|
| **C** | baked table, **per-index** fallback above it | baked table, same fallback | 93.1 MB |
| **C+D** | baked table + fallback | **host hints, checked `t·h==1`** | 50.1 MB |
| **M** | **host leaves, one pruned multi-path vs a baked root** | host hints, checked | **7.15 MB** |

All three remove the pre-2026-09-13 cliff. A hinted *inverse* is safe because it is checkable in one
multiply; a hinted *twist* is not — it selects which code is verified against — so M commits to it.

## Results — instructions and segments

| sigs | C | C+D | M |
|---|---|---|---|
| 32 | 67,522,115 (7) | 71,751,216 (7) | 86,424,091 (9) |
| 64 | 72,609,342 (7) | 77,202,200 (8) | 94,918,684 (10) |
| 128 | 78,610,104 (8) | 83,679,713 (9) | 105,508,607 (11) |
| 256 | 89,116,798 (9) | 94,215,932 (10) | 115,900,587 (13) |
| 512 | 109,809,029 (11) | **103,825,476 (10)** | 126,322,449 (15) |
| 1024 | 131,247,719 (12) | **114,167,254 (12)** | 137,796,349 (17) |

**D's sign flips with batch size.** Below the table it costs ~5 M (a checked hint lookup against a
34-instruction table read); above it, it *saves* — C still pays a full field inversion at each
derived level, D replaces that with one multiply. At 1024 that is **−17 M**. D is the batch-inversion
fix in disguise, not only an image optimisation.

## Results — cells (`--mode meter`, `Total cost`)

| sigs | C | C+D | M |
|---|---|---|---|
| 64 | 2,881,612,925 | 3,055,564,714 | 3,741,012,958 |
| 256 | 3,523,137,015 | 3,716,554,512 | 4,556,589,626 |
| 512 | 4,389,673,359 | **4,097,893,138** | 4,963,285,519 |
| 1024 | 5,285,770,738 | **4,508,280,555** | 5,411,723,571 |

## Results — `prove app`

Same box, same session, `--exe` throughout.

| sigs | C | C+D | M |
|---|---|---|---|
| 32 | 25 s / 5,537,700 | **23 s** / 5,600,566 | 24 s / 6,940,754 |
| 128 | 26 s / 6,583,535 | **24 s** / 6,878,914 | 30 s / 8,812,288 |
| 256 | 28 s / 7,040,789 | 29 s / 7,860,935 | **27 s** / 9,472,897 |
| 1024 | 38 s / 10,723,886 | **29 s** / 9,391,572 | 30 s / 11,742,583 |

**At 1024, C+D is 9 s faster than C** — D's 17 M instruction saving showing up as wall-clock, with a
smaller app proof too (9.39 vs 10.72 MB).

Wall-clock repeats within ±2 s (C+D here vs the image-cost run: 23/24/29/29 against 23/24/29/27), so
differences under ~2 s are noise. The 9 s at 1024 is not.

## Blob size

| sigs | C | C+D | M |
|---|---|---|---|
| 32 | 1,457,927 | 1,597,327 (+9.6 %) | 2,011,843 (+38 %) |
| 256 | 2,215,003 | 2,391,203 (+8.0 %) | 3,051,479 (+38 %) |
| 1024 | 2,630,259 | 2,846,299 (+8.2 %) | 3,814,991 (+45 %) |

Aggregation erases this later — the root proof is a constant 534,636 B — so it is a prover-side
cost, not a permanent artefact size.

## Verdict

**C+D.** Best instructions and cells at every size above 256, fastest `prove app` at the largest
batch by 9 s, a smaller app proof than either alternative at 1024, half the image of C, and only 8 %
more blob.

**M is a specialised tool.** It delivers what nothing else can — a **7.15 MB** image, 13x smaller,
and it removes `MAX_LOG_COLS` as a concept rather than softening it. But it costs +23 % instructions,
+5 segments, +38–45 % blob and +20–28 % app proof, to save ~43 MB of image that
`2026-09-13-recursion-image-size-cost` measured at **~0.04 s/MB — about 1.7 s**. Take it only when
image size is a hard constraint (a deployment limit, a slow artefact channel), not for speed.

**Had the documented 0.25 s/MB stood, M would have looked like the winner.** That 43 MB would have
read as ~11 s of savings instead of ~1.7 s. Measuring the image cost directly is what inverted the
conclusion.

## Comparability

- All three columns: same box (vultr4), same commit, same signature counts, same invocation form,
  same session. Blobs regenerated in the current format for all three.
- **not comparable to:** any `prove app` before 2026-09-11 (build path included, ~18 s), or to the
  first M attempt (crashed on the `g0` and roots-array bugs; no numbers taken).

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| Merkle-committing the twist table (M) | works; 13x smaller image | +23 % instructions, +5 segments, +38–45 % blob, for ~1.7 s of image saving |
| growing `MAX_LOG_COLS` instead of C | not attempted | image is `2^k`; C gets the same reach for 46 k instructions and 6 KB |
| sizing `TWIST_ROOTS` by `MAX_LOG_COLS` | guest panic at 1024 | 1024 needs `log_cols` 18; roots must be sized independently |
| building M's payload from recorded indices alone | guest panic at index 0 | `g0` is read by the code *constructor*, never through `diag_at`, so the recorder cannot see it |

## Notes

**M moves cost from the image into build time**: `build.rs` samples ~8.4 M elements and builds 21
trees single-threaded — **811 s**. Parallelising would fix it; not done.

**M's pruning advantage erodes with scale**: 6.5x at tree height 16, 4.6x at height 20. Leaves grow
~1.5x across the range while the tree grows 16x, so queries share fewer ancestors — the wrong
direction for a strategy meant to scale.

**C's overhead below the table is ~46 k instructions (0.06 %)** — the bounds check. That is the whole
price of removing the cliff, and it is why C should be adopted regardless of what happens to D or M.
