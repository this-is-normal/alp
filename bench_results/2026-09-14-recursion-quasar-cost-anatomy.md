# Anatomy of the Quasar recursion guest: **cryptography is 37 %**, and 60 % of the guest is deserialisation plus opening machinery

<!-- layer: recursion -->
<!-- target: Quasar recursion guest, full phase x operation decomposition -->

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | where the guest's instructions actually go, decomposed on **two independent axes** |
| **parameters** | `num_bits=253`, S=256 and 4096, strategy C+D, `MAX_LOG_COLS=16`, `BF_INV_RATE=4`, `POW_BITS=16` |
| build features | guest `hint-inverses` + one `prof-stop-*` or `ablate-*` per cut · host `quasar,bb-hash,parallel,count-ops` + the same `ablate-*` |
| commands | `cargo openvm run --exe … --mode segment` (guest, per cut) · `count_verify_bf 253 <S>` (host, per cut) |

## Method

Two measurements, deliberately independent:

1. **Phase cost — measured in-guest.** Nine `ablate-*` cuts plus two new `prof-stop-*` cuts make the
   verifier return early; the delta between adjacent cuts is that phase's real instruction cost.
2. **Cryptographic content — exact counts, host-side.** The same cuts applied to
   `count_verify_bf`, so `count-ops` reports field and hash operations **per phase** rather than
   globally.

Dividing one by the other gives each phase's *accounting ratio*: how much of it is cryptography.

**New for this run:** `prof-stop-deser` / `prof-stop-config` on the guest (with `black_box`, or the
compiler elides the work being measured), and six `ablate-*` passthroughs on `eddsa-example` which
previously exposed only three.

## Results — S=4096, guest total 135,042,075

| phase | measured | cryptography | accounted | non-crypto | % of guest |
|---|---|---|---|---|---|
| **deserialisation** | 33,911,700 | 0 | **0 %** | 33,911,700 | 25.1 % |
| config construction | 67,650 | 0 | — | 67,650 | 0.1 % |
| zerocheck | 2,367,164 | 1,918,052 | 81 % | 449,112 | 1.8 % |
| **sampled-column check** | 31,720,546 | 27,423,375 | **86 %** | 4,297,171 | 23.5 % |
| QA-encoding sumchecks | 4,315,966 | 4,315,966 | 106 %* | 0 | 3.2 % |
| claim reduction | 766,944 | 766,944 | 172 %* | 0 | 0.6 % |
| **both Basefold openings** | 61,958,420 | 15,018,347 | **24 %** | 46,940,073 | 45.9 % |
| closing AIR check | 878,316 | ~500 | ~0 % | 877,821 | 0.7 % |
| **CRYPTOGRAPHY** | | **49,443,179** | | | **36.6 %** |
| **NON-CRYPTOGRAPHIC** | | | | **86,543,527** | **64.1 %** |

\* over 100 %: the op-count x microbench model is an **upper bound**. Capped at measured in the
totals. Small phases overshoot most, having less to average over.

Phase costs sum to 0.7 % over the measured total; per-phase op counts sum to 50,273,115 against the
50,273,110 measured independently — the two axes close.

## The headline

**Cryptography is 37 % of the guest, and it sits where there is least leverage.**

| | measured | crypto | grows with batch? |
|---|---|---|---|
| column check | 31.7 M | **86 %** | **no** — 30.6 -> 31.7 M across 16x |
| openings | 62.0 M | **24 %** | **yes** — 70 % of all growth |

The phase that is genuinely cryptographic is fixed cost. The phase that scales is three-quarters
machinery: transcript bookkeeping, index arithmetic and `Vec` handling around the fold-chain walk.

**The non-cryptographic 64 % has two blocks:**

| target | instructions | share of non-crypto |
|---|---|---|
| **opening machinery** | 46.9 M | **54 %** |
| **deserialisation** | 33.9 M | **39 %** |
| everything else | 5.7 M | 7 % |

Together **80.9 M — 60 % of the guest, containing essentially no cryptography.**

## Scaling: the guest becomes *less* cryptographic as the batch grows

| | S=256 | S=4096 |
|---|---|---|
| cryptography | 45.9 % | **36.6 %** |
| non-cryptographic | 54.9 % | **64.1 %** |
| openings, accounted | 32 % | **24 %** |

Because the crypto-heavy column check is fixed while the machinery-heavy openings scale.

## Two measured constants

**Deserialisation costs 10.1 instructions per blob byte** — 10.12 at S=256, 10.14 at S=4096, across
a 1.4x blob range. Against an input floor of 3.24 instr/byte (`read_vec` + touching every byte),
an addressable proof format has a ceiling of ~6.9 instr/byte, i.e. **~23 M at 4096, ~17 % of the
guest** — more if regions can be left untouched entirely.

**Config construction is 68 k instructions**, flat across a 16x batch range — 0.05 % of the guest,
from 309 M (56 %) before the parameters were baked and the resolved shape shipped.

## Comparability

- **comparable to:** `2026-09-13-recursion-quasar-cd-scaling` — same box, commit, build, blobs; its
  totals are the denominators here.
- **not comparable to:** the 2026-09-11 phase attribution, which predates strategies C and D.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| the Basefold-internal cuts (`ablate-after-sumcheck`, `ablate-after-merkle`) | **discarded** | `bf_sumcheck` lands 12.5 M *below* the `reduction` cut that must precede it — the documented signature of a cut removing code that exists elsewhere. Their difference (51.4 M) is suggestive but unverifiable |
| the tool's Merkle row (`8,071 paths x 9,593`) | **discarded** | 77.4 M exceeds the whole verify and double-counts permutations already itemised; the 9,593 constant was measured at depth 14 with 136 sparse queries and does not transfer to dense ones |
| single-axis attribution | insufficient | the operation axis alone said "crypto is a minority"; the phase axis alone said "openings dominate". Both are true and neither is actionable without the other |

## Notes

**The documented 3.4x paired-opening artifact did not reproduce.** This run measures the two
openings at 31.5 M and 30.5 M — a ratio of **0.97x**. The pair (62.0 M) is now corroborated by
near-equal halves rather than resting on the counter argument alone. The trap remains real; it is
simply not present in the current code shape.

**`before_close` sits 878 k below `both_opens`** — not an ordering violation. They cut in different
layers: the Quasar cut returns from the PCS and multi-stark still runs its closing AIR check, which
the multi-stark cut skips. The gap matches `full - before_close` to within 1.3 k.

**What this implies for the work queue.** Addressable proofs target the second-largest
non-cryptographic block (33.9 M, ceiling ~17 %). The largest — **opening machinery at 46.9 M**, and
the only one that grows — has no proposal against it yet. It is a question about the verifier's data
structures and loop shape, not about Quasar.
