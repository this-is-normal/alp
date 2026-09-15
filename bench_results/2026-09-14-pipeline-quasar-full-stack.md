# Full stack S=32..4096: 4096 signatures -> a 534,636 B proof at 8.6 ms/sig — and the **native prover becomes the bottleneck above ~2,700**

<!-- layer: pipeline (prover + recursion + aggregation) -->
<!-- target: Quasar, strategy C+D, all three layers -->

## Target

| | |
|---|---|
| layer | `pipeline` — native prover, `prove app`, `prove stark` |
| subject | Quasar with **strategy C+D** (hybrid table lookup + host-hinted checked inverses) |
| **parameters** | `num_bits=253`, S=32..4096, `BF_INV_RATE=4`, `POW_BITS=16`, `MERKLE_CAP_HEIGHT=0`, `MAX_LOG_COLS=16`, `--segment-max-memory 26 GiB` |
| build features | guest `hint-inverses` · host `quasar,bb-hash,parallel(,record-queries)` |
| input artefacts | `cd2_*.bin`, current blob format (trailing `merkle_twists` element) |
| `.vmexe` | 50,115,207 B |
| **box load at start** | **3.07** — not idle; see Notes |

## 1) Native Plonky3 prover (Quasar)

| sigs | trace | proof B | B/sig | setup ms | prove ms | ms/sig | verify ms |
|---|---|---|---|---|---|---|---|
| 32 | 8192x45 | 1,428,701 | 44,646.9 | 81.0 | 954.5 | 29.83 | 110.9 |
| 64 | 16384x45 | 1,610,534 | 25,164.6 | 157.0 | 1,543.0 | 24.11 | 107.1 |
| 128 | 32768x45 | 1,822,813 | 14,240.7 | 236.6 | 2,528.9 | 19.76 | 106.3 |
| 256 | 65536x45 | 2,196,130 | 8,578.6 | 255.9 | 3,792.3 | 14.81 | 124.2 |
| 512 | 131072x45 | 2,411,232 | 4,709.4 | 397.4 | 6,848.6 | 13.38 | 133.6 |
| 1024 | 262144x45 | 2,670,273 | 2,607.7 | 679.1 | 12,963.3 | 12.66 | 147.6 |
| 2048 | 524288x45 | 2,925,839 | 1,428.6 | 1,160.9 | 24,794.2 | **12.11** | 160.1 |
| 4096 | 1048576x45 | **3,202,008** | **781.7** | 1,927.1 | 52,977.0 | 12.93 | **172.2** |

Prove time per signature bottoms out at **12.11 ms** around 2048. **Verify runs 110.9 -> 172.2 ms
across a 128x batch range** — the polylog verifier behaving as designed.

## 2) `prove app` and 3) `prove stark`

3 repeats for app, 2 for stark; medians with spread.

| sigs | instr M | seg | cells M | app s | ± | app proof B | stark s | ± | **agg s** | root B |
|---|---|---|---|---|---|---|---|---|---|---|
| 32 | 71.8 | 7 | 2,843 | 22.8 | 0.6 | 5,600,566 | 23.6 | 0.2 | 0.79 | 534,636 |
| 64 | 77.2 | 8 | 3,056 | 26.0 | 0.9 | 6,368,538 | 26.8 | 0.8 | 0.78 | 534,636 |
| 128 | 83.7 | 9 | 3,309 | 24.4 | 0.6 | 6,878,914 | 24.4 | 2.0 | **−0.01** | 534,636 |
| 256 | 94.2 | 10 | 3,717 | 28.1 | 0.6 | 7,860,935 | 29.6 | 0.4 | 1.50 | 534,636 |
| 512 | 103.9 | 10 | 4,099 | 27.7 | 1.1 | 8,359,249 | 29.3 | 0.2 | 1.58 | 534,636 |
| 1024 | 114.2 | 12 | 4,509 | 27.4 | 0.2 | 9,391,572 | 29.3 | 0.7 | 1.93 | 534,636 |
| 2048 | 124.6 | 13 | 4,920 | 31.7 | 1.4 | 10,298,456 | 33.7 | 0.1 | 1.99 | 534,636 |
| 4096 | 135.0 | 14 | 5,336 | 33.6 | 1.2 | 11,347,093 | 35.4 | 0.2 | 1.82 | 534,636 |

**The root proof is exactly 534,636 B at all eight sizes.**

### Aggregation does NOT scale with segments — measurement debt paid

`prove stark` costs **0.8-2.0 s** on top of `prove app`. Across a doubling of segments (7 -> 14) it
goes 0.79 -> 1.82 s, and the 128-signature point measured **−0.01 s**, i.e. inside the noise.

Treat aggregation as **~1.5 s flat**. The previous estimate ("+2 s", from two points at 1 s
resolution) was numerically close but rested on no evidence about scaling.

## End-to-end per signature

| sigs | native ms/sig | app ms/sig | stark ms/sig | root B/sig |
|---|---|---|---|---|
| 32 | 29.83 | 711.9 | 736.6 | 16,707 |
| 256 | 14.81 | 109.8 | 115.7 | 2,088 |
| 1024 | 12.66 | 26.8 | 28.7 | 522 |
| **4096** | **12.93** | **8.20** | **8.64** | **130.5** |

**4096 ed25519 signatures compress to 534,636 bytes at 8.6 ms and 130 bytes per signature**, both
still falling with batch size.

## THE HEADLINE: the bottleneck has moved

| sigs | native | stark | total | native share | bottleneck |
|---|---|---|---|---|---|
| 32 | 1.0 s | 23.6 s | 24.6 s | 4 % | recursion |
| 256 | 3.8 s | 29.6 s | 33.4 s | 11 % | recursion |
| 1024 | 13.0 s | 29.3 s | 42.3 s | 31 % | recursion |
| 2048 | 24.8 s | 33.7 s | 58.5 s | 42 % | recursion |
| **4096** | **53.0 s** | 35.4 s | 88.4 s | **60 %** | **NATIVE** |

The two layers scale differently and the gap widens without bound:

- **native: linear in S** — it proves a trace that grows with the batch (~12.5 ms/sig);
- **recursion: logarithmic in S** — it verifies *one proof*, +1.7 s per doubling.

**Crossover is ~2,700 signatures.** At 8192 the native prover would be ~102 s against ~37 s.

**Consequence.** Every optimisation in this project to date targets the recursion guest, which above
2,700 signatures is the *smaller and shrinking* half. Below that it remains correct. Two mitigations
keep this from being as bad as it looks:

1. **Different hardware.** Native proving is CPU (rayon, 30.72-core cgroup quota); `prove app` is
   GPU. A pipeline over successive batches overlaps them, so throughput is `max(native, recursion)`
   ~ 53 s rather than the 88 s sum.
2. **Headroom.** 53 s on a 30-core quota with 256 vCPUs visible, and rayon measured 8.1x on ~23
   cores. The GPU path already saturates one 5090.

## Comparability

- **comparable to:** `2026-09-13-recursion-quasar-cd-scaling` — same box, commit, guest and blobs;
  app figures agree within the repeat spread.
- **not comparable to:** any wall-clock before 2026-09-11 (build path included, ~18 s), or runs on a
  loaded box (load average was ~22 during some earlier sweeps).

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| treating aggregation as segment-proportional | not supported | 7 -> 14 segments moved it 0.79 -> 1.82 s, with one point inside the noise |

## Notes

**Box load was 3.07**, not idle. App spreads run 0.2-2.0 s and stark 0.1-2.0 s; differences below
~2 s are not claimable. The negative aggregation figure at 128 signatures is the clearest
illustration.

**The app-proof curve is nearly flat** — 22.8 s at 32 signatures against 33.6 s at 4096, a 128x
batch for 1.5x the time. That flatness is what collapses the per-signature cost, and it argues for
running the largest batch memory allows.
