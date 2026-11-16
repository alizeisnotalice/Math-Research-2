# Completed n=12 GPU inverse-fibre batch

This is a finite one-peak SMC diagnostic for canonical inverse fibres of distance 2 and 3, with `n=12`, `P=65536`, `R=128`, `tau=145`, and `L=16`. It estimates `A2_distance_d/A1`; it is neither an exact n=12 count nor an estimate of the original two-peak conditional distribution.

`inverse_fibre_gpu_r020.swift` is the actual production source, SHA-256 `e56997103ba08ea5507b667499f3d059349c6c4b96c00cf1599484bdf78f2e56`. The raw JSONL, sidecar metadata, production summary, device telemetry, and independent root reaggregation are included. The production binary and traces are intentionally omitted.

```sh
python3 summarize_n12_production.py
python3 verify_n12_archive.py
```

The first command regenerates the production summary from raw JSONL. The second cross-checks its common estimates with `n12-inverse-root-summary.json`, validates the recorded production source/raw hashes, and checks the completed telemetry exit record.

`P65536_SAFETY.md` gives the exact scope of the uint indexing audit. It certifies this production P value, not every P allowed by the old general entry guard.

The copied production metadata and frozen manifest preserve one historical absolute `baseline_source` string as evidence of their original run context. Reproduction does not depend on that path: the actual e569 source and all raw inputs used by the local scripts are in this directory.
