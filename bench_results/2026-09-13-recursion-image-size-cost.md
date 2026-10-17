# `.vmexe` size costs ~0.04 s/MB of `prove app`, not 0.25 — the documented figure is ~6x too high

<!-- layer: recursion -->
<!-- target: OpenVM program-image size vs prove-app wall-clock -->

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | does baked static data in the program image cost proving time, and how much |
| **parameters** | `num_bits=253`, S=32/128/256/1024, `MAX_LOG_COLS=16`, `--segment-max-memory 26 GiB` |
| build features | guest: `hint-inverses` (small) vs none (large) · host blob `quasar,bb-hash,parallel,record-queries` |
| command | `cargo openvm prove app --exe <vmexe> --input <json> --segment-max-memory $((26<<30))` |
| input artefacts | `hint_quasar_*.bin` — **the same blobs for both arms** |

## Why this is a clean isolation

Both builds read blobs that carry inverse hints, and the guest's
`hinted_twist_inverse` prefers hints whenever they are present. So **both arms
execute the same instructions and emit the same proof**; they differ only in
whether the now-unused 43 MB `TINV_WORDS` table is baked into the image.

That the app proofs come out **byte-identical in every pair** confirms it.

## Results

| sigs | SMALL 50,089,128 B | LARGE 93,080,472 B | delta | app proof (both) |
|---|---|---|---|---|
| 32 | 23 s | 25 s | +2 s | 5,600,566 |
| 128 | 24 s | 26 s | +2 s | 6,878,914 |
| 256 | 29 s | 31 s | +2 s | 7,983,494 |
| 1024 | 27 s | 28 s | +1 s | 9,391,572 |

**42,991,344 bytes of image costs 1–2 s: ~0.023–0.047 s/MB, call it 0.04.**

## What this corrects

`alp-openvm/MANIFEST.md` and `OPTIMIZATIONS.md` both record that *"a 92 MB
`.vmexe` costs ~16 s more per `prove app` than a 27.7 MB one on identical
work"* — 0.25 s/MB. At that rate the 43 MB here should have cost ~10.7 s. It
costs ~2.

The likely cause is the same contamination found on 2026-09-11: wall-clock
measured with `cargo openvm prove app` **including the build/keygen path**,
which is worth ~18 s on its own. Two builds of different sizes measured that way
would differ by however much that path differed, not by the image cost. This run
uses `--exe` throughout, so no build path is in either arm.

**Consequences, all of which were being reasoned about with the wrong number:**

- Deleting `TINV_WORDS` (strategy D) is worth **~2 s**, not ~10 s. Its real value
  is elsewhere: it removes the field inversion at derived levels, worth
  **−17 M instructions at 1024 signatures**.
- The argument against growing `MAX_LOG_COLS` to 18 (~179 MB image) weakens
  considerably: ~3.4 s, not ~21 s. It is still the worst lever — it is charged at
  every batch size and strategy C gets the same benefit for 46 k instructions —
  but the case against it is now about instruction cost, not proving time.
- Strategy M's headline benefit — removing ~86 MB of image — is worth **~3.4 s**,
  against a blob 38 % larger. That is a much weaker case than the image figure
  suggested.

## Comparability

- **comparable to:** each other only. Both arms: same box, same session, same
  blobs, same `--exe` invocation, back to back.
- **not comparable to:** any wall-clock measured with the build path in the loop,
  including every `prove app` figure in this repo before 2026-09-11.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| trusting the recorded 0.25 s/MB | wrong by ~6x | it was a modelled attribution from two runs that differed in more than image size |

## Notes

**Baked static data is nearly free at proving time.** The original claim that it
is "free at execution and free at proving up to ~2 MB" was directionally right
and too conservative on the ceiling: 43 MB costs ~2 s.

**This does not make image size irrelevant** — it is charged at every batch size,
so a 179 MB image still taxes the 32-signature case that would never benefit. But
image size should no longer be the *primary* argument in any of these trade-offs.
