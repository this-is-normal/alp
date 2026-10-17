---
name: bench-record
description: Record a benchmark result in bench_results/ with full provenance. Use when the user has just run or asked to record a benchmark of a PCS, the Plonky3 prover, an OpenVM recursion guest, or OpenVM aggregation — or says "record this bench", "write up these numbers", "add to bench_results".
---

# Recording a benchmark

Write one `.md` narrative and one `.json` sidecar into `bench_results/`, then regenerate the index.
The format and its rationale are in `bench_results/README.md`; the shape is `bench_results/TEMPLATE.md`.

**This skill records results. It does not run benchmarks.** The runs are too varied — host sweeps,
guest instruction counts, GPU proves — to hide behind one command.

## Steps

1. **Collect provenance mechanically. Never from recall or memory.**

   ```bash
   ./bench_results/bench_record.sh            # measured locally
   ./bench_results/bench_record.sh <host>     # measured on a remote box
   ```

   Paste the output verbatim into the Provenance section. If the run was on a remote box, pass that
   host — the CPU, GPU and CUDA must describe **where the numbers came from**, not where you are.

2. **Fill the Target table completely.** Every parameter that could change the result, not just the
   interesting one. If you don't know one, find it in the source rather than omitting it — a
   constant buried in a build script has cost 2.7× here without appearing in any table.

3. **Put hardware-independent numbers first**: instructions, segments, cells, field-op counts,
   proof bytes. Wall-clock second, always naming the box.

4. **Label every number `measured` or `(modelled)`.** Anything from a cost model is an estimate.
   Op-count × microbench-cost is an **upper bound** — it has overshot by 6×, because a
   `black_box`ed loop forces spills real code elides. Never present a modelled number as a result.

5. **Fill Comparability explicitly**, both directions. A run is not comparable to another with
   different parameters, a different box (for wall-clock), or a commit that touched the hot path.
   If the numbers supersede an earlier record, say so **in that record too**.

6. **Fill "Tried and rejected", even if empty.** Negative results are results, and without a home
   they get re-attempted.

7. Name the files `YYYY-MM-DD-<layer>-<target>-<variant>.{md,json}`, then:

   ```bash
   python3 bench_results/make_index.py
   ```

## Refuse to finish if

- provenance was typed by hand rather than produced by the script
- the **uncommitted** count is non-zero and the record does not say so — "commit abc123" is a lie
  when the tree has 25 uncommitted files, which has happened here
- target parameters are missing or partial
- measured and modelled numbers are mixed without labels
- a comparison is drawn against a run whose provenance differs, without the caveat spelled out

## Sanity checks worth making before writing

- Does the input artefact (proof blob, trace) come from **this** commit? A stale blob has looked
  exactly like a logic bug before.
- Were any `ablate-*` or `count-ops` features on? They change codegen and instrument the code —
  say so, and never quote such a run as a headline.
- If wall-clock is being compared across boxes, is a hardware-independent metric quoted alongside?
