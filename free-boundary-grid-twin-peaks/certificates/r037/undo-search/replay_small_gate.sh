#!/bin/bash
# Public small gate for r046. It expects r036 as a sibling certificate.
set -euo pipefail

base=$(cd "$(dirname "$0")" && pwd)
set +u
baseline_src=$BASELINE_SOURCE
steps=$STEPS
set -u
if [ -z "$baseline_src" ]; then baseline_src="$base/../../r036/finite-candidates/psi_search.cpp"; fi
if [ -z "$steps" ]; then steps=10000; fi
expected_baseline=98fcb4848a23230d533cff96f68ed4fd8162ddbb96223877bae1211e1719debc
expected_undo=4f4be34e7ea8066dfe00817f1dbc60d039bd45566d52f602afc27c451d5200e6
out="${REPLAY_OUT:?set REPLAY_OUT to a fresh data-volume directory}"
build="$out/build"
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang

test -f "$baseline_src" || {
  echo "missing frozen r040 source: $baseline_src (set BASELINE_SOURCE to override)" >&2
  exit 2
}
test "$(shasum -a 256 "$baseline_src" | awk '{print $1}')" = "$expected_baseline"
test "$(shasum -a 256 "$base/psi_search_undo.cpp" | awk '{print $1}')" = "$expected_undo"
mkdir -p "$out" "$build"

c++ -O2 -std=c++20 -pthread "$baseline_src" -o "$build/psi_search_r040"
c++ -O2 -std=c++20 -pthread "$base/psi_search_undo.cpp" -o "$build/psi_search_undo"

run_case() {
  local n=$1
  local root=$2
  local seed=$3
  local init=$4
  local tag="n$n-r$root-s$seed-$init"
  "$build/psi_search_r040" "$n" "$root" "$seed" "$steps" "$out/r040-$tag.txt" 1 "$init"
  "$build/psi_search_undo" "$n" "$root" "$seed" "$steps" "$out/undo-$tag.txt" 1 "$init"
  python3 "$base/verify_psi.py" "$out/undo-$tag.txt" >"$out/undo-$tag.verify"
  cmp "$out/r040-$tag.txt" "$out/undo-$tag.txt"
  printf 'exit_status=0\ncomparison=byte_identical\n' >"$out/undo-$tag.receipt"
}

run_case 3 4 1 bfs
run_case 3 0 2 ring
run_case 4 10 3 bfs
run_case 4 0 1 row
printf 'exit_status=0\ncases=4\nproposals_per_case=%s\naudit_interval=1\n' "$steps" >"$out/ALL_SMALL_GATES.receipt"
echo "PASS: 4 cases, audit=1, output=$out"
