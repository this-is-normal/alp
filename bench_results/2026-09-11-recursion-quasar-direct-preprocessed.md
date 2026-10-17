# Quasar recursion guest, direct preprocessed evaluation: 62.3 M instr @ 16 sigs (10.4x since 2026-09-09); a cost cliff above 256 signatures

<!-- layer: recursion -->
<!-- target: Quasar recursion guest (OpenVM) -->

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | Quasar PCS, verified inside an OpenVM guest |
| **parameters** | `num_bits=253`, signatures 16..1024, `c=2`, `BF_INV_RATE=4`, `POW_BITS=16`, `MERKLE_CAP_HEIGHT=0`, **`MAX_LOG_COLS=16`** |
| build features | `quasar,bb-hash,parallel` (host) · `rv32i,rv32m,io,modular,poseidon2` chips (guest) |
| command | `cargo openvm run --exe <vmexe> --input <json> --mode segment` · `cargo openvm prove app --exe <vmexe> --input <json> --segment-max-memory $((26<<30))` |
| input artefacts | `gen_proof 253 <n>` regenerated at this commit, 2026-09-11 13:21-13:39; blob sizes in the results table |
| what changed vs the previous record | the preprocessed trace is **no longer committed**: `MultiStarkConfig::evaluates_preprocessed_directly()` lets the verifier evaluate its MLE in closed form, because the EdDSA preprocessed trace is one `h`-row pattern tiled `2^k` times |

## Provenance

| field | value |
|---|---|
| date | 2026-09-11 |
| tree | `alp2` @ `70cb077d` (`quazar2`) |
| **uncommitted** | **true (19 files)** |
| host | vultr4 |
| CPU | AMD EPYC 7B12 64-Core Processor — 256 visible, quota 30.72 cores |
| memory | 503 GB |
| GPU | NVIDIA GeForce RTX 5090, 32607 MiB, 580.159.03 |
| rustc | 1.98.1 |
| cargo-openvm | cargo-openvm v2.0.2 (59a69b8) [cuda] |

The 19 uncommitted files **are** the change under test. This record is not reproducible from
`70cb077d` alone.

## Results

Hardware-independent first. All variants reveal `[1, 0, ...]`.

| sigs | blob bytes | instructions | segments | code table | `prove app` (vultr4) | app proof bytes | s/sig |
|---|---|---|---|---|---|---|---|
| 16 | 1,305,070 | 62,343,966 | 6 | baked | 21 s | 5,202,291 | 1.31 |
| 32 | 1,457,918 | 67,476,735 | 7 | baked | 24 s | 5,537,700 | 0.75 |
| 64 | 1,619,622 | 72,561,795 | 7 | baked | 24 s | 5,781,529 | 0.375 |
| 128 | 1,806,390 | 78,563,621 | 8 | baked | 27 s | 6,706,670 | 0.211 |
| 256 | 2,214,994 | 89,069,686 | 9 | baked | 28 s | 7,028,421 | 0.109 |
| 512 | 2,404,410 | **345,930,649** | **31** | **lazy fallback** | 89 s | 23,081,846 | 0.174 |
| 1024 | 2,630,250 | **367,371,064** | **33** | **lazy fallback** | 89 s | 24,295,107 | **0.087** |

All figures **measured**. The `prove app` column uses the `--exe` form throughout, which skips the
build/keygen path; see Notes. Every row in it was taken with the identical invocation, so the
column is internally comparable.

Below the cliff the wall-clock is nearly flat — **21 s at 16 signatures to 28 s at 256**, a 16x
batch for a 1.3x cost — which is the whole point of recursion as a compression play.

### The cliff above 256 signatures

`MAX_LOG_COLS=16` sizes the baked twist+inverse table. The chooser gives `N2 = 2^16` at 256 (fits),
`2^17` at 512 and `2^18` at 1024 (both overflow), so those two fall back to
`RandomFoldableCode::new_verifier`, which derives twists at the queried positions and pays a field
inversion for each.

The fallback is a **step, not a slope**: 512 -> 1024 adds only +6% (21.4 M instructions) for a 2x
batch, because doubling adds exactly one code level. Backing that out:

| | |
|---|---|
| marginal cost per code level | 21.4 M instructions |
| inversion cost (microbench) | 14,879 instructions |
| implied queried positions per level | ~1,440 |
| levels at `2^18` | 18 |
| implied total inversion cost | ~386 M — i.e. essentially the whole excess |

So ~95% of the fallback's cost is one operation. Montgomery batch inversion turns `q` inversions
into one inversion plus ~`3q` multiplications at 66 instructions each — ~200 per position instead
of 14,879. *(modelled)* 512 and 1024 would land near ~100 M, restoring the ~10 M-per-doubling trend
that 16..256 shows.

Raising `MAX_LOG_COLS` to 18 is **not** the fix: the table is `2^k`, so the `.vmexe` goes from
93 MB to ~370 MB, and 93 MB already costs ~16 s per `prove app`.

### Native prover, same change

| | before | after |
|---|---|---|
| setup @ 1024 sigs | 9,736 ms | 528 ms |
| prove @ 1024 sigs | 14,070 ms | 12,386 ms |
| verify @ 1024 sigs | 261 ms | 147 ms |
| proof @ 1024 sigs | 4.86 MB | 2.67 MB |

## Comparability

- **comparable to:** `2026-09-10-recursion-quasar-flat` — same box, same parameters, same chips;
  the only difference is the preprocessed commitment removal. 16 sigs 87.1 M -> 62.3 M.
- **not comparable to:** that record's **wall-clock** column — it used `cargo openvm prove app`
  with the build path in the loop; this one uses `--exe`. See Notes.
- **not comparable to:** `2026-09-09-recursion-basefold` — BaseFold predates this and several other
  shared-crate fixes that would lift it too.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| raising `MAX_LOG_COLS` to cover `2^18` | not attempted — rejected on size | table is `2^k`; 93 MB `.vmexe` -> ~370 MB, and 93 MB already costs ~16 s per `prove app` |

## Notes

**The wall-clock column is not comparable to earlier records.** These runs pass `--exe` and skip
the build/keygen path; 32/64/128 previously measured 42/45/43 s *with* it. The ~16 s that
`MANIFEST.md` attributes to loading a 93 MB `.vmexe` is about the whole difference, which is why
256's 28 s must not be read as faster than 128's 43 s. The matched re-measure settles it: 32 goes
42 s -> 24 s and 128 goes 43 s -> 27 s on identical work, so the overhead was ~18 s and 256 is in
fact a hair *slower* than 128 (28 s vs 27 s), as it should be.

**512 is the worst operating point in the sweep.** It pays the entire fixed cliff and amortises it
over half the signatures of 1024. If deploying above 256 today, go to 1024.

**Two tooling failures preceded this run, both silent-looking.** `--input` takes OpenVM's
`{"input": ["0x01<hex>"]}` wrapper, not the raw blob, and the leading `01` is a type tag. The first
repair used `xxd`, which is **not installed on vultr4**, so every JSON was 19 bytes and the CLI
rejected the empty hex with a message about the type tag. The encoder is now validated by
regenerating a known-good `bf_input_16.json` and comparing byte-for-byte.

**A stale `.app.proof` reported a plausible size through six failed invocations.** The driver now
`rm -f`s the proof before each attempt so a failure reports `MISSING`. Third time this class of
error has cost time here.
