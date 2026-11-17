# Replay and audit

This package needs only a C++20 compiler, `shasum`, POSIX `sh`, and Python 3
standard library.  It has no third-party Python dependencies.

From the package root, audit the full shipped finite evidence without compiling
or rescanning any graph.  First pre-create a writable directory *outside* this
frozen certificate; the auditor's optional JSON output belongs there:

```sh
OUT=/path/to/precreated/n9-audit-output
test -d "$OUT"
python3 scripts/audit_n9_certificate.py --root . --json-out "$OUT/n9-audit-result.json"
```

That command reads the 326 small checkpoint JSON files, recomputes their
SHA-256 values, checks their source SHA and zero integer relation, reconstructs
the strict intervals, checks the 326 exit records, and validates both gates.
Omit `--json-out` to print the result only; neither form writes under the
certificate root.

For a fresh minimum executable replay, choose a disposable writable directory
outside this certificate and a shared writable compiler cache:

```sh
export CLANG_MODULE_CACHE_PATH=/path/to/shared/clang-cache
sh scripts/min_replay.sh /path/to/disposable/n9-min-replay
python3 scripts/audit_min_replay.py \
  --work /path/to/disposable/n9-min-replay \
  --source-sha be18e9f8c593648c4e2121305f400db5199ed2fcdde46c4522ddbae27b15b5c4
```

The minimum replay compiles the source, runs the N=7 compatibility gate, runs
the three-graph N=9 DP/direct-`9!` gate, and scans exactly 16 left-1 host edge
subsets.  It intentionally does not rerun the complete 1,327,360-host batch.
The actual package-preparation replay passed with Apple clang 21.0.0 and Python
3.15.0; see `evidence/min-replay-receipt.json`.  A regenerated executable may
have a different binary hash on a different compiler/platform.

For a new full run, use the source orchestration only after supplying a
deliberate output directory and `CLANG_MODULE_CACHE_PATH`; it is outside this
certificate's replay scope.  Do not run it into `evidence/`, whose files are
frozen evidence.
