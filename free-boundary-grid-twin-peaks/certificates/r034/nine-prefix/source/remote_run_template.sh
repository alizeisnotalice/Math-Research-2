#!/bin/sh
# Template for the single authorized SSH agent. DO NOT run this on the local Mac.
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
OUT="$ROOT/results"
SRC="$ROOT/bipartite_prefix_dp.cpp"
BIN="$ROOT/bipartite_prefix_dp"
SHA=$(shasum -a 256 "$SRC" | awk '{print $1}')
# Fresh-only launch. Existing receipts are evidence: do not reuse this output
# directory. Use resume_aggregate.py to prepare a distinct suffix shard.
if [ -e "$OUT" ]; then
  echo "refusing to overwrite existing results directory: $OUT" >&2
  exit 2
fi
mkdir -p "$OUT"

# Compile and gate remotely. Any exit 3 means preserve the JSON counterexample
# and do not launch or continue the affected full shard.
c++ -O3 -std=c++20 -Wall -Wextra -Werror -pedantic "$SRC" -o "$BIN"
"$BIN" --mode gate-n7-r45 --out-dir "$OUT" --source-sha "$SHA"
"$BIN" --mode gate-n9-direct --out-dir "$OUT" --source-sha "$SHA"
python3 "$ROOT/make_shard_plan.py"

# The coordinator launches at most fourteen workers. Redirect each worker log;
# inspect every shard JSON. A status other than completed is never aggregated.
python3 - "$ROOT/shard-plan.json" "$BIN" "$OUT" "$SHA" <<'PY'
import json, subprocess, sys
plan, binary, out, sha = sys.argv[1:]
procs=[]
for s in json.load(open(plan))["shards"]:
    cmd=["nice","-n","10",binary,"--mode","shard-r27","--left-size",str(s["left_size"]),"--start",str(s["start"]),"--stop",str(s["stop"]),"--shard-id",s["shard_id"],"--out-dir",out,"--source-sha",sha,"--checkpoint-every","1000"]
    log=open(f"{out}/{s['shard_id']}.log","wb")
    procs.append((s,log,subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)))
bad=[]
for s,log,p in procs:
    rc=p.wait(); log.close()
    if rc: bad.append((s["shard_id"],rc))
if bad: raise SystemExit(f"non-completed shards (preserved for inspection): {bad}")
PY
