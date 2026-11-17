# Replaying the completed evidence

Use existing Python3 and a C++20 compiler. Choose a fresh output directory
outside the certificate tree and an existing shared compiler cache:

```sh
TMPDIR=/Volumes/data/CodexTemp/T2 \
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang \
PYTHONDONTWRITEBYTECODE=1 \
python3 -B free-boundary-grid-twin-peaks/certificates/r034/replay.py \
  --out /Volumes/data/CodexTemp/T2/a022-r034-new-replay
```

This runs four Python checks, authenticates all 326 stored N9 checkpoints,
and compiles the original C++ source to replay the N7 gate, the three-graph
N9 direct-permutation cross-check, and a 16-graph shard. It does not repeat the
complete 1,327,360-graph production run. The detailed evidence-only audit and
coverage scope are in nine-prefix/REPLAY.md.

Scripts that write beside themselves are copied into external work directories.
Mathematical JSON is compared exactly; only the interpreter-version field in
the singleton-image JSON is excluded. Ring computations use exact Fraction
inequalities internally; their saved decimals are descriptive displays.

parent-replay-receipt.json records the actual final-tree execution. Every child
was waited for and exited zero. This replay does not stop or authenticate the
separate in-progress production jobs, which are outside the round34 results.
