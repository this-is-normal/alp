#!/bin/bash
# Reproduce the batched EdDSA comparison in QUASAR_BENCHMARK_REPORT.md §9.
#
# Sweeps every commitment backend over signature counts, printing one line per
# (backend, count). Parse with the snippet at the bottom of this file to
# regenerate `bench_data/eddsa_batch_parallel.json`.
#
#   ./run_batch_sweep.sh                      # parallel, BabyBear hash (the report's §9)
#   HASH= PARALLEL= ./run_batch_sweep.sh      # serial, ed25519 hash (§5)
#   COUNTS=1,4,16 ./run_batch_sweep.sh        # smaller sweep
#   PROFILE=1 ./run_batch_sweep.sh            # also print the prover span tree
#
# Run from the workspace root. Set CARGO_TARGET_DIR to a local (non-sshfs) path:
# the mount's linker cannot link binaries in place.
set -u
: "${COUNTS:=1,4,16,64,256}"
: "${NUM_BITS:=253}"
: "${HASH:=bb-hash}"
: "${PARALLEL:=parallel}"
: "${CARGO:=cargo}"

for backend in "" ligero basefold quasar; do
  feats=$(IFS=,; echo "${backend:+$backend,}${HASH:+$HASH,}${PARALLEL}" | sed 's/,$//;s/^,//')
  echo "### backend=${backend:-brakedown} features=${feats:-none}"
  $CARGO run -q -p eddsa-example --release --bin measure \
    ${feats:+--features "$feats"} -- "$NUM_BITS" "$COUNTS" 2>&1 \
    | grep -aE "^(Brakedown|Ligero|BaseFold|Quasar|error)"
done

# --- Parse the output into JSON -------------------------------------------------
# ./run_batch_sweep.sh | tee sweep.txt
# python3 - <<'PY'
# import re, json
# rows = []
# for line in open('sweep.txt'):
#     m = re.match(r'(\w+)\s+sigs\s+(\d+)\s+trace=\s*(\d+)x(\d+)\s+proof\s+(\d+) B '
#                  r'\(\s*[\d.]+ B/sig\)\s+setup\s+([\d.]+) ms\s+prove\s+([\d.]+) ms '
#                  r'\(\s*([\d.]+) ms/sig\)\s+verify\s+([\d.]+) ms', line.strip())
#     if m:
#         rows.append(dict(backend=m[1], signatures=int(m[2]), trace_rows=int(m[3]),
#                          trace_cols=int(m[4]), proof_bytes=int(m[5]), setup_ms=float(m[6]),
#                          prove_ms=float(m[7]), prove_ms_per_sig=float(m[8]),
#                          verify_ms=float(m[9])))
# json.dump(dict(rows=rows), open('bench_data/eddsa_batch_parallel.json', 'w'), indent=1)
# PY
