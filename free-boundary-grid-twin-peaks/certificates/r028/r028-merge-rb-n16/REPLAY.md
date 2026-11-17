# Portable finite receipt replay

The production raw is included uncompressed. The trace JSON is losslessly
gzip-compressed because its uncompressed size exceeds the per-file archive
limit. `UNCOMPRESSED_SHA256.txt` records the required decompressed SHA-256.

Run from this directory, setting `TRACE_WORK` to a writable external working
directory.

```sh
gzip -dc trace/n16-p17-r1-l16.jsonl.trace.json.gz > "$TRACE_WORK/n16-trace.json"
shasum -a 256 "$TRACE_WORK/n16-trace.json"
python3 -B source/validate_trace_payload.py \
  trace/n16-p17-r1-l16.jsonl "$TRACE_WORK/n16-trace.json" \
  source/tau2_merge_rb_calibration_r026.swift > verify-payload.json
cmp reports/n16-p17-r1-l16.payload-audit.json verify-payload.json
python3 -B source/verify_merge_rb_trace.py \
  trace/n16-p17-r1-l16.jsonl "$TRACE_WORK/n16-trace.json" > verify-rank.json
cmp reports/n16-p17-r1-l16.rank-replay.json verify-rank.json
python3 -B source/aggregate_production_r026.py \
  raw/n16-p262144-r64-l16.jsonl source/tau2_merge_rb_calibration_r026.swift > verify-production.json
cmp reports/n16-production-audit.json verify-production.json
```

The r028 validator repair also checks that both trace identifiers, rather
than only a chained inequality expression, are zero.  Its four-case replay
uses this included n=16 trace: `0/0` must pass and `1/1`, `0/1`, and `1/0`
must each fail.

```sh
python3 -B trace-id-audit/test_trace_ids.py \
  trace-id-audit/validate_trace_payload_r028.py \
  trace/n16-p17-r1-l16.jsonl "$TRACE_WORK/n16-trace.json" \
  source/tau2_merge_rb_calibration_r026.swift \
  "$TRACE_WORK/id-validation-test-replay.json"
python3 -B trace-id-audit/verify_test_results.py \
  "$TRACE_WORK/id-validation-test-replay.json"
```

These checks validate a finite trace path, raw-file input invariants, and
derived arithmetic. The production audit explicitly does not certify SMC
accuracy, mixing, cone bounds, or asymptotics.
