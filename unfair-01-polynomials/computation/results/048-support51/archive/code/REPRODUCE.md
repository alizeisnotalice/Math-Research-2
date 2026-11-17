# Semantic total-51 cache coverage audit

This is a bounded cap-six audit for **total degree 51 only**.  It establishes
that the retrieved cache has one passed, hash-bound audited result for every
independently generated endpoint-normalized support `I` in each degree pair
`m + n = 51`, with `1 <= m <= 25`.  It does not replay every remote compressed
chunk ledger and therefore is not a second semantic replay of the classification.

The immutable copied inputs are in `input/`; their expected SHA-256 values are
hard-coded in `audit51.py` and emitted in `report51.json`.  The cache input is
`cache/t51/` (68,405 JSON records).  The cache relative content list digest
is also emitted in the report.

To reproduce after transferring the same remote files into this directory:

```sh
nice -n 10 /Users/zzh/.local/bin/python3 \
  /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/semantic51_audit47/audit51.py
```

The script checks byte bindings, independently regenerates the exact 68,405
support/index table, validates complete and duplicate-free event coverage,
matches every cache record to its event and frozen binding, and checks
`progress51.json` plus `exit51.json` for successful joined completion.

For the separate one-core raw-corpus audit (already run in this evidence
directory), execute:

```sh
nice -n 10 /Users/zzh/.local/bin/python3 \
  /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/semantic51_audit47/raw51_bytes_audit.py
```

It reads every local `round042-macbook-support52/run/t51/*.jsonl.gz` file,
checks the compressed and decompressed SHA-256 values against the event, checks
newline row count against the event's canonical count, and rechecks the
matching cache binding.  It deliberately does not perform a second row-wise
semantic classification.
