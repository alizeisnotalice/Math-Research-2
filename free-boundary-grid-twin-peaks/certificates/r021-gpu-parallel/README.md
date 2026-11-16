# r021 optimized terminal inverse-fibre kernel

This minimal archive freezes `r021_optimized.swift` at SHA-256
`6aa1f829dbe97c141781c7440c0f99ef62fa966cbe2c665fab7432a48b2ec9f8`.
It supports `2<=n<=16`, ordered prefixes of stride 256, and the rank-parallel
terminal d=2,3 inverse-fibre calculation.  The particle limit is checked
before allocation so `particle*256+vertex` remains a Metal `uint` address.

The source preserves two historical output-schema fields deliberately:
the raw JSONL calls itself r021 and records terminal-integration base source
SHA `1f593425af5fe53f98cc3b01109d6498519ad2d555c2ba14332d715c07f26b48`,
while the trace JSON retains the frozen r020 trace-schema identifier and its
SHA.  `verify_trace.py` and `verify_inverse.py` check the latter schema and
reconstruct the graph process independently; they do not silently equate the
two source identities.

`n3-final-l16.jsonl` plus its 176 KiB trace is an actual r021 run.  The
independent inverse result is `n3-final-l16.jsonl.inverse-host.json`.
`original-final-validation.md` records the source's final n=3 agreement and
the pre-allocation invalid-particle rejection.

## Reproduce from this directory

```sh
TASK_TMP=/Volumes/data/CodexTemp/T2/a022-r021-archive-replay
mkdir -p "$TASK_TMP"
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules \
  r021_optimized.swift -o "$TASK_TMP/r021"
"$TASK_TMP/r021" 3 17 1 10 10 "$TASK_TMP/n3.jsonl" 1f2e3d4c strat 16 trace
python3 verify_trace.py "$TASK_TMP/n3.jsonl.trace.json"
python3 verify_inverse.py "$TASK_TMP/n3.jsonl.trace.json"
```

The programs write all replay output below `TASK_TMP`; no archive file is
overwritten.  This is an implementation audit, not a conditional-probability
estimator.
