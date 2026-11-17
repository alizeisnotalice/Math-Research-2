# Portable replay

Run from this directory, with `TRACE_WORK` set to a writable external cache.

```sh
gzip -t raw/n24-p1048576-r32-l16.jsonl.gz
gzip -dc raw/n24-p1048576-r32-l16.jsonl.gz > "$TRACE_WORK/n24.jsonl"
shasum -a 256 "$TRACE_WORK/n24.jsonl"
python3 -B source/frozen_r023_aggregate_tau2_distance.py "$TRACE_WORK/n24.jsonl" \
  --sidecar raw/aggregate.sidecar.json --source source/tau2_distance_terminal_r023.swift \
  --epsilon 1 > verify-epsilon-1.json
python3 -B source/frozen_r023_aggregate_tau2_distance.py "$TRACE_WORK/n24.jsonl" \
  --sidecar raw/aggregate.sidecar.json --source source/tau2_distance_terminal_r023.swift \
  --epsilon 1/2 > verify-epsilon-1_2.json
python3 -B source/frozen_r023_aggregate_tau2_distance.py "$TRACE_WORK/n24.jsonl" \
  --sidecar raw/aggregate.sidecar.json --source source/tau2_distance_terminal_r023.swift \
  --epsilon 1/4 > verify-epsilon-1_4.json
cmp replay-epsilon-1.json verify-epsilon-1.json
cmp replay-epsilon-1_2.json verify-epsilon-1_2.json
cmp replay-epsilon-1_4.json verify-epsilon-1_4.json
python3 -B audit_raw_and_cone.py "$TRACE_WORK/n24.jsonl" raw/aggregate.sidecar.json \
  source/tau2_distance_terminal_r023.swift verify-epsilon-1_2.json verify-audit.json
cmp raw-and-cone-audit.json verify-audit.json
python3 -B certify_cone_discrepancy.py verify-epsilon-1_2.json verify-discrepancy.json
cmp cone-discrepancy-certificate.json verify-discrepancy.json
```

The uncompressed SHA must equal the `raw_sha256` entry in
`raw/RAW_GZIP_SHA256.txt`.  The replays establish finite-file integrity and
deterministic arithmetic only; see `README.md` for the scope limitation.
