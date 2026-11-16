# r021 n<=32 extension

`n32_extension.swift` has SHA-256
`fc45ad181f6de297001c835a77d9be4620b944169269ed3cdb0f9c9a9a5a8629`.
It is self-contained and changes prefix/rank stride to 1024 while retaining
the exact 32-word occupancy and boundary state.  It requires `P<=2^22` before
allocation, so every `particle*1024+vertex` index is within Metal `uint`.

`n32-one-final-p1-l2.jsonl` and its 9.1 MiB trace are the smallest full n=32
one-peak audit.  They demonstrate a non-dead Double underflow: finite `logZ`,
zero `Z` and `Zmean_*`, but nonzero `mean_*`/`scaled_mean_*`.  The inverse
verifier checks this representation explicitly.

`p17_outputs/` contains the four small raw P=17 outputs and resource records.
Their 52--152 MiB traces are deliberately omitted.  `p17_trace_manifest.json`
records their SHA-256 values, commands and validation status; the separately
run root replay receipt is `p17_root_replay_receipt.json`.

`run_negative_trace_tests.py` changes only disposable copies of the included
P=1 trace in `TASK_TMP` (or `/Volumes/data/CodexTemp/T2` by default), and
removes them after each run.  Its saved `negative-tests.json` proves rejection
of a modified high word, truncated prefix and changed inverse group.

## Reproduce from this directory

```sh
TASK_TMP=/Volumes/data/CodexTemp/T2/a022-n32-archive-replay
mkdir -p "$TASK_TMP"
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules \
  n32_extension.swift -o "$TASK_TMP/n32"
"$TASK_TMP/n32" 32 1 1 1025 1025 "$TASK_TMP/n32.jsonl" 3f4e5d6c strat 2 trace
python3 verify_trace.py "$TASK_TMP/n32.jsonl.trace.json"
python3 verify_inverse.py "$TASK_TMP/n32.jsonl.trace.json"
python3 run_negative_trace_tests.py
```

The last command uses the archived P=1 trace and creates only disposable
negative copies under TASK_TMP.  None of these results estimates a
conditional peak probability.
