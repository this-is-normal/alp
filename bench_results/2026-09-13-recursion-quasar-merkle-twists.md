# Strategy M — Merkle-committed twist table: image 93 MB -> 7.15 MB, but +23 % instructions and +5 segments

<!-- layer: recursion -->
<!-- target: Quasar recursion guest, twist table replaced by a baked Merkle root -->

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | replacing the baked twist table with a baked **root**; the host ships queried leaves with one pruned multi-path |
| **parameters** | `num_bits=253`, S=32..1024, `BF_INV_RATE=4`, `POW_BITS=16`, `MERKLE_CAP_HEIGHT=0`, `MAX_ROOT_LOG_COLS=20` |
| build features | guest `merkle-twists,hint-inverses` · host `quasar,bb-hash,parallel,record-queries` + `MERKLE_TWISTS=1` |
| command | `cargo openvm run --exe <vmexe> --input <json> --mode segment` / `--mode meter` |
| input artefacts | `merkle_*.bin`, regenerated 2026-09-13 after the `g0` fix |

## The design

The twists are public and fixed, so committing to them is sound: the guest checks the host's leaves
against a root in its **own program image**, covered by the exe commitment exactly as the table was.
Unlike an inverse, a twist is *not* checkable in one multiply — it selects which code is being
verified against — so a commitment is the only way to accept one from the host.

One root per `log_cols` (21 of them, 672 B) rather than one big tree, so path length and
verification cost match the shape. One tree serves both the main and preprocessed codes: same seed,
prefix index layout.

## Results

All reveal `[1, 0, …]`. **`.vmexe` 7,153,272 B** (was 93,042,506 — **13x smaller**).

| sigs | instructions | segments | cells | blob bytes |
|---|---|---|---|---|
| 32 | 86,424,091 | 9 | — | 2,011,843 |
| 64 | 94,918,684 | 10 | 3,741,012,958 | 2,297,675 |
| 128 | 105,508,607 | 11 | 4,154,175,721 | 2,649,179 |
| 256 | 115,900,587 | 13 | 4,556,589,626 | 3,051,479 |
| 512 | 126,322,449 | 15 | 4,963,285,519 | 3,407,615 |
| 1024 | 137,796,349 | 17 | 5,411,723,571 | 3,814,991 |

### Tree shape, size and pruning

| sigs | leaves | tree | **height** | pruned siblings | path bytes | unpruned siblings | pruning gain |
|---|---|---|---|---|---|---|---|
| 32 | 3,489 | 2^16 | 16 | 8,592 | 274,944 | 55,824 | 6.5x |
| 64 | 3,893 | 2^17 | 17 | 11,461 | 366,752 | 66,181 | 5.8x |
| 128 | 4,419 | 2^18 | 18 | 15,294 | 489,408 | 79,542 | 5.2x |
| 256 | 4,409 | 2^18 | 18 | 15,122 | 483,904 | 79,362 | 5.2x |
| 512 | 4,929 | 2^19 | 19 | 19,032 | 609,024 | 93,651 | 4.9x |
| 1024 | 5,405 | 2^20 | 20 | 23,515 | 752,480 | 108,100 | 4.6x |

**Height tracks `log_cols`, not signature count** — 128 and 256 share height 18 because they share a
code shape, which is also why their leaf counts are nearly equal.

**Pruning's advantage erodes as the tree grows**: 6.5x at height 16 down to 4.6x at height 20. Leaves
grow ~1.5x across the range while the tree grows 16x, so queries get sparser and share fewer
ancestors. That is the wrong direction for a strategy meant to scale.

## Comparison at the same sizes

| | image | blob @256 | instr @256 | seg @256 | instr @1024 | seg @1024 |
|---|---|---|---|---|---|---|
| C | 93.0 MB | 2,214,994 | 89,116,798 | 9 | 131,247,719 | 12 |
| **C+D** | 50.1 MB | 2,391,202 | 94,215,932 | 10 | **114,167,254** | **12** |
| **M** | **7.15 MB** | 3,051,479 | 115,900,587 | 13 | 137,796,349 | 17 |

M against C+D: **−43 MB image, +21.7 M instructions (+23 %), +3 segments, +0.66 MB blob** at 256;
at 1024, +23.6 M (+21 %) and +5 segments.

Since `2026-09-13-recursion-image-size-cost` measured image at **~0.04 s/MB**, that 43 MB is worth
**~1.7 s** — and segments drive `prove app` directly. See that record for the `prove app` outcome.

## Comparability

- **comparable to:** `2026-09-13-recursion-quasar-hybrid-table` (C) and the C+D sweep — same box,
  same commit, same signature counts, same `--mode segment`/`meter` invocation.
- **not comparable to:** the first M attempt on 2026-09-13, which used blobs generated before the
  `g0` fix and a roots array sized by `MAX_LOG_COLS`. Both crashed; no numbers were taken.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| sizing `TWIST_ROOTS` by `MAX_LOG_COLS` | guest panic, `len is 17 but the index is 18` | 1024 signatures need `log_cols` 18; roots are 32 B each and must be sized independently of the table |
| building the payload from recorded indices alone | guest panic, "must ship a twist for every queried index: 0" | `g0` (`0..inv_rate`) is read by the code's **constructor**, not through `diag_at`, so the recorder never sees it |

## Notes

**Strategy M moves cost from the image into build time.** `build.rs` now samples ~8.4 M code
elements and builds 21 trees, single-threaded: **811 s**. The image shrinks 13x and the generated
file goes from 23 MB to 2,521 B, but a clean build costs 13.5 minutes. Parallelising the sampling
would fix it; it has not been done.

**The instruction cost is the pruned verification**, ~21.7 M at 256 for 15,122 sibling compressions
plus 4,409 leaf hashes. That is far cheaper per path than the 9,593 instructions a single pruned
query-path costs in isolation — dense queries share almost everything — but it is not free, and it
grows with tree height while the table's cost does not.

**The blob grows 38–45 %**, which aggregation later erases (the root proof is a constant 534,636 B),
so this is a prover-side cost, not a permanent artefact-size one.
