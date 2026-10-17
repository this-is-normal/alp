#!/usr/bin/env bash
# Print the provenance block for a benchmark record.
#
# Exists because the mistakes worth preventing are omissions, not bad
# arithmetic: a run compared against a different scalar size, a different box,
# or a tree with uncommitted changes. Collect it mechanically, never by recall.
#
#   ./bench_results/bench_record.sh [--json] [remote-host]
#
# With a host argument the hardware section is collected over ssh, for runs
# measured on a remote box.
set -uo pipefail

FORMAT=md
[ "${1:-}" = "--json" ] && { FORMAT=json; shift; }
HOST="${1:-}"

# Remote shells are non-interactive, so cargo/CUDA are not on PATH by default.
run() {
  local cmd="source ~/.cargo/env 2>/dev/null; export PATH=\$PATH:/usr/local/cuda/bin; $1"
  if [ -n "$HOST" ]; then ssh -o BatchMode=yes "$HOST" "$cmd" 2>/dev/null; else eval "$cmd" 2>/dev/null; fi
}

repo_root=$(git rev-parse --show-toplevel 2>/dev/null || echo .)
commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)
branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)
dirty_n=$(git status --porcelain 2>/dev/null | wc -l | tr -d ' ')
dirty=$([ "$dirty_n" = "0" ] && echo false || echo "true ($dirty_n files)")

host_name=${HOST:-$(hostname)}
cpu=$(run 'lscpu | sed -n "s/^Model name: *//p"' | head -1)
cores=$(run 'nproc')
quota=$(run 'cat /sys/fs/cgroup/cpu.max 2>/dev/null' | awk '{ if ($1=="max") print "none"; else printf "%.2f cores", $1/$2 }')
mem=$(run 'free -g | awk "/^Mem/{print \$2\" GB\"}"')
gpu=$(run 'nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader' | head -1)
cuda=$(run 'nvcc --version 2>/dev/null | tail -1' | sed -n 's/.*release \([0-9.]*\).*/\1/p')
rustc=$(run 'rustc --version' | awk '{print $2}')
openvm=$(run 'cargo openvm --version' | head -1)

if [ "$FORMAT" = json ]; then
  printf '{\n'
  printf '  "date": "%s",\n' "$(date -u +%Y-%m-%d)"
  printf '  "commit": "%s", "branch": "%s", "dirty": "%s",\n' "$commit" "$branch" "$dirty"
  printf '  "host": "%s", "cpu": "%s", "cores": "%s", "cpu_quota": "%s",\n' "$host_name" "$cpu" "$cores" "$quota"
  printf '  "memory": "%s", "gpu": "%s", "cuda": "%s",\n' "$mem" "$gpu" "$cuda"
  printf '  "rustc": "%s", "cargo_openvm": "%s"\n' "$rustc" "$openvm"
  printf '}\n'
else
  echo "| field | value |"
  echo "|---|---|"
  echo "| date | $(date -u +%Y-%m-%d) |"
  echo "| tree | \`$(basename "$repo_root")\` @ \`$commit\` (\`$branch\`) |"
  echo "| **uncommitted** | **$dirty** |"
  echo "| host | $host_name |"
  echo "| CPU | $cpu — $cores visible, quota $quota |"
  echo "| memory | $mem |"
  [ -n "$gpu" ] && echo "| GPU | $gpu |"
  [ -n "$cuda" ] && echo "| CUDA | $cuda |"
  echo "| rustc | $rustc |"
  [ -n "$openvm" ] && echo "| cargo-openvm | $openvm |"
fi
