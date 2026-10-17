# Strategy C+D to 4096 signatures: **8.3 ms/signature**, and why `prove app` barely moves across a 32x batch range

<!-- layer: recursion -->
<!-- target: Quasar recursion guest, strategy C+D, scaling to 4096 -->

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | strategy C+D (hybrid table lookup + host-hinted checked inverses) across S=128..4096 |
| **parameters** | `num_bits=253`, `BF_INV_RATE=4`, `POW_BITS=16`, `MERKLE_CAP_HEIGHT=0`, `MAX_LOG_COLS=16`, `--segment-max-memory 26 GiB` |
| build features | guest `hint-inverses` · host `quasar,bb-hash,parallel,record-queries` |
| command | `cargo openvm run --exe … --mode segment` / `--mode meter`; `cargo openvm prove app --exe … --segment-max-memory $((26<<30))` |
| input artefacts | `cd2_*.bin`; **2048 and 4096 generated for the first time** (193 s and 385 s natively) |
| `.vmexe` | 50,115,207 B |

**Three `prove app` repeats per size.** An earlier sweep showed a non-monotonic result (more
segments proving *faster*), implying >=3 s of noise. This box is a shared vast.ai container and its
load average was ~22 during that sweep; it was **1.47** during this one.

## Results

| sigs | instructions | seg | cells | `prove app` reps | mean | spread | app proof | ms/sig |
|---|---|---|---|---|---|---|---|---|
| 128 | 83,703,303 | 9 | 3,309,227,735 | 25, 24, 25 | 24.7 s | 1 s | 6,878,914 | 192.7 |
| 256 | 94,239,522 | 10 | 3,717,296,426 | 30, 28, 29 | 29.0 s | 2 s | 7,860,935 | 113.3 |
| 512 | 103,852,440 | 10 | 4,098,745,912 | 27, 27, 28 | 27.3 s | 1 s | 8,359,249 | 53.4 |
| 1024 | 114,197,592 | 12 | 4,509,244,189 | 28, 28, 27 | 27.7 s | 1 s | 9,391,572 | 27.0 |
| 2048 | 124,555,808 | 13 | 4,919,842,858 | 32, 32, 32 | 32.0 s | **0 s** | 10,298,456 | 15.6 |
| 4096 | 135,042,075 | 14 | 5,336,055,212 | 34, 33, 35 | 34.0 s | 2 s | 11,347,093 | **8.3** |

**Repeat spreads are 0–2 s**, against the >=3 s and an inverted ordering seen on the loaded box. The
earlier wall-clock was contaminated by co-tenants; these numbers are usable, those were not.

## Why the curve is almost flat

Three effects compound, and only the first is about the scheme:

**1. The guest verifies one proof, not N signatures.** Work grows **logarithmically** in the batch:

| step | instructions added |
|---|---|
| 128 -> 256 | +10.5 M |
| 256 -> 512 | +9.6 M |
| 512 -> 1024 | +10.3 M |
| 1024 -> 2048 | +10.4 M |
| 2048 -> 4096 | +10.5 M |

**+10.4 M per doubling on an ~84 M base** — dead linear in `log(S)`. Across 32x the signatures the
guest does **1.61x** the work. This is the whole point of proof recursion, and it is the dominant
term.

**2. Segments grow slower still** — 9 -> 14 across that 32x, tracking cells (1.61x).

**3. The cost is quantised.** 256 and 512 have the **same 10 segments** and the same ~28 s, though
512 executes 10 % more instructions. `prove app` steps with segment count, and the steps are coarse:
9 -> ~25 s, 10–12 -> ~28 s, 13 -> 32 s, 14 -> 34 s.

So "128 and 1024 take the same time" was an overstatement of a real effect: they differ (24.7 vs
27.7 s), but by far less than 8x, because the guest's work is logarithmic and the prover's cost is
a coarse step function of segments.

## Comparability

- **comparable to:** `2026-09-13-recursion-strategy-comparison`'s C+D column — same box, commit,
  build and invocation. Small deltas against it (128: 83,703,303 here vs 83,679,713 there, +23,590)
  are the blob's extra trailing `Option` element, which that sweep's blobs predated.
- **not comparable to:** wall-clock from the earlier sweeps on the loaded box.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| single-sample `prove app` timings | produced a non-monotonic ordering (M at 256 faster than at 128 with more segments) | on a shared box the noise exceeded the signal; repeats were needed |

## Notes

**4096 is the first real stress of strategy C's fallback.** The chooser puts `log_cols` around 20,
so ~4 levels sit above the baked `MAX_LOG_COLS=16` table and are derived per index. The instruction
curve stays linear in `log(S)` right through it — the per-index fallback degrades gracefully, which
is exactly what the all-or-nothing version did not do.

**One mild inversion survives: 256 (29.0 s mean) is slower than 512 (27.3 s)**, despite the same
segment count and fewer instructions. The ranges barely separate ([28,30] vs [27,28]). Not
explained; possibly uneven segment fill. Recorded rather than smoothed over.

**Per-signature cost falls 23x across the range**, 192.7 ms at 128 to **8.3 ms at 4096**, and is
still falling — the fixed cost has not yet been amortised out at 4096.
