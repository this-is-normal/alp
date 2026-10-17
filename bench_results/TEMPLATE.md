# <one line: what was measured and the headline number>

<!-- layer: pcs | prover | recursion | aggregation -->
<!-- target: e.g. Quasar recursion guest -->

## Target

| | |
|---|---|
| layer | `recursion` |
| subject | Quasar PCS, verified inside an OpenVM guest |
| **parameters** | `num_bits=253`, signatures 1..128, `c=2`, `BF_INV_RATE=4`, `POW_BITS=16`, `MAX_LOG_COLS=16` |
| build features | `quasar,bb-hash,parallel` (host) · `modular,poseidon2` chips (guest) |
| command | `cargo openvm prove app --input … --segment-max-memory $((26<<30))` |
| input artefacts | blob from `gen_proof 253 <n>`, regenerated at this commit |

Every parameter that can change the result belongs here, not just the interesting one. A constant
buried in a build script has cost 2.7× before.

## Provenance

<!-- paste ./bench_results/bench_record.sh [host] verbatim -->

## Results

Hardware-independent first.

| variant | instructions | segments | proof bytes | wall-clock (box) |
|---|---|---|---|---|
| | | | | |

All figures **measured** unless marked *(modelled)*.

## Comparability

- **comparable to:** `<run>` — same parameters, same box, only `<X>` differs
- **not comparable to:** `<run>` — `<reason>`

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| | | |

Empty is fine; absent is not.

## Notes

Anything a reader needs to reproduce or to avoid misreading the numbers.
