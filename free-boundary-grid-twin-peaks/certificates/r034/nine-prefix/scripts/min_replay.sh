#!/bin/sh
# Compile once and replay two gates plus one 16-graph shard; never the full scan.
set -eu
if [ "$#" -ne 1 ]; then
  echo "usage: $0 WORK_DIR" >&2
  exit 64
fi
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
WORK=$1
case "$WORK" in "") echo "WORK_DIR must be nonempty" >&2; exit 64;; esac
mkdir -p "$WORK/gates" "$WORK/results"
BIN="$WORK/bipartite_prefix_dp"
SRC="$ROOT/source/bipartite_prefix_dp.cpp"
SHA=$(shasum -a 256 "$SRC" | awk '{print $1}')
: "${CLANG_MODULE_CACHE_PATH:?set this to a shared writable compiler cache}"
export CLANG_MODULE_CACHE_PATH
c++ -O3 -std=c++20 -Wall -Wextra -Werror -pedantic "$SRC" -o "$BIN"
"$BIN" --mode gate-n7-r45 --out-dir "$WORK/gates" --source-sha "$SHA"
"$BIN" --mode gate-n9-direct --out-dir "$WORK/gates" --source-sha "$SHA"
"$BIN" --mode shard-r27 --left-size 1 --start 0 --stop 16 --shard-id min-l1-0-16 --out-dir "$WORK/results" --source-sha "$SHA" --checkpoint-every 8
python3 "$ROOT/scripts/audit_min_replay.py" --work "$WORK" --source-sha "$SHA"
