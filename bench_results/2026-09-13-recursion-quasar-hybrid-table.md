# Strategy C — per-index table fallback: the cliff above 256 signatures is gone (3.15x at 512, 2.80x at 1024)

<!-- layer: recursion -->
<!-- target: Quasar recursion guest, hybrid baked-table lookup -->

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | replacing the all-or-nothing baked-table lookup with a **per-index** fallback |
| **parameters** | `num_bits=253`, S=32..1024, `BF_INV_RATE=4`, `POW_BITS=16`, `MERKLE_CAP_HEIGHT=0`, `MAX_LOG_COLS=16` (**unchanged**) |
| build features | none (default) · host blob `quasar,bb-hash,parallel` |
| command | `cargo openvm run --exe <vmexe> --input <json> --mode segment` and `--mode meter` |
| input artefacts | the same `eddsa_quasar_*.bin` blobs as the 2026-09-11 records, unchanged |

## The change

Guest-side only; **no alp2 change, no blob change, no table growth**. `baked_twist` and
`baked_twist_inverse` gained a per-index fallback:

```rust
if idx < twists::TWIST_WORDS.len() { table[idx] } else { sample_nonzero_at(seed, idx) }
```

The index layout is a prefix (`0..inv_rate` for `g0`, then `level_base(inv_rate, level) + j`), so a
table sized for `MAX_LOG_COLS` covers levels `1..=MAX_LOG_COLS` **completely** whatever `log_cols`
the proof uses. The old code abandoned the entire table when one index fell outside it, so 1024
signatures derived all **18** levels; this derives **2**.

## Results

All reveal `[1, 0, …]`. `.vmexe` 93,042,506 B (was 93,036,034 — the +6 KB is the fallback code).

| sigs | instructions before | **instructions (C)** | delta | segments before | **segments (C)** | cells (`Total cost`) |
|---|---|---|---|---|---|---|
| 32 | 67,476,735 | 67,522,115 | +45,380 | 7 | 7 | 2,682,461,761 |
| 64 | 72,561,795 | 72,609,342 | +47,547 | 7 | 7 | 2,881,612,925 |
| 128 | 78,563,621 | 78,610,104 | +46,483 | 8 | 8 | 3,116,215,550 |
| 256 | 89,069,686 | 89,116,798 | +47,112 | 9 | 9 | 3,523,137,015 |
| **512** | 345,930,649 | **109,809,029** | **−236,121,620** | 31 | **11** | 4,389,673,359 |
| **1024** | 367,371,064 | **131,247,719** | **−236,123,345** | 33 | **12** | 5,285,770,738 |

**512: 3.15x. 1024: 2.80x.** Segment counts, which drive `prove app`, fell 31->11 and 33->12.

**The cost below the table is the bounds check alone: ~+46 k instructions, 0.06 %.** That is the
entire price of removing the cliff.

The curve is now smooth — 67.5, 72.6, 78.6, 89.1, 109.8, 131.2 M — with no discontinuity. The
steeper slope above 256 is the two derived levels, at ~14.8 M each, exactly as measured.

### The cost model predicted this

`2026-09-11-recursion-quasar-lazy-cost-attribution` measured **14.76 M instructions per derived
level** and projected C at ~113 M (512) and ~140 M (1024).

| | projected | measured | error |
|---|---|---|---|
| 512 | ~113 M | 109,809,029 | −2.8 % |
| 1024 | ~140 M | 131,247,719 | −6.3 % |

Both came in slightly **better** than modelled. Two independent measurements now support the
per-level figure.

## Comparability

- **comparable to:** `2026-09-11-recursion-quasar-direct-preprocessed` — same box, same blobs, same
  parameters; the only difference is this change. Its rows are the "before" column above.
- **not comparable to:** anything measured before 2026-09-11 — the preprocessed-commitment removal
  sits in between.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| growing `MAX_LOG_COLS` to 17 or 18 instead | not attempted | table is `2^k`; measured image relation is `7.03 MB + 5.12 x data`, so 2^18 is ~179 MB of `.vmexe` — charged at **every** batch size, to help only the largest. C costs 46 k instructions and 6 KB |

## Notes

**`MAX_LOG_COLS` is no longer a cliff, but it is still a slope.** Shapes above it cost ~14.8 M
instructions per excess level. The knob now trades image size against a gentle per-level cost
instead of falling off a 3x edge, so mis-sizing it is survivable.

**This does not reduce the image.** Both tables are still baked; `.vmexe` is still ~93 MB. That is
what strategy D (hinted inverses, deleting `TINV_WORDS`) addresses, separately.

**The cells column** (`cargo openvm run --mode meter`, reported as `Total cost`) is recorded here
as the hardware-independent proving-cost proxy, for comparison against the Merkle-tree strategy.
