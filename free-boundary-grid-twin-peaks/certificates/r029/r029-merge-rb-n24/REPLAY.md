# Portable finite receipt replay

Run from this directory, setting `A022_REPLAY_DIR` to a writable cache.  The
two compressed inputs are lossless; their compressed and uncompressed hashes
are recorded in `INPUT_SHA256.txt`.

```sh
gzip -dc raw/n24-p262144-r32-l16.jsonl.gz > "$A022_REPLAY_DIR/raw.jsonl"
gzip -dc trace/n24-p17-r1-l16-trace.json.gz > "$A022_REPLAY_DIR/trace.json"
python3 -B source/aggregate_production_r026.py "$A022_REPLAY_DIR/raw.jsonl" source/tau2_merge_rb_calibration_r026.swift > "$A022_REPLAY_DIR/audit.json"
cmp reports/n24-production-audit.json "$A022_REPLAY_DIR/audit.json"
python3 -B source/validate_trace_payload.py trace/n24-p17-r1-l16.jsonl "$A022_REPLAY_DIR/trace.json" source/tau2_merge_rb_calibration_r026.swift > "$A022_REPLAY_DIR/payload.json"
cmp reports/n24-p17-r1-l16.payload-audit.json "$A022_REPLAY_DIR/payload.json"
python3 -B source/verify_merge_rb_trace.py trace/n24-p17-r1-l16.jsonl "$A022_REPLAY_DIR/trace.json" > "$A022_REPLAY_DIR/rank.json"
cmp reports/n24-p17-r1-l16.rank-replay.json "$A022_REPLAY_DIR/rank.json"
```

This verifies one finite implementation receipt, not mixing, SMC accuracy,
the conditional law, or an asymptotic claim.
