# Completed r022 tau=2 calibration

This subpackage is the completed calibration snapshot from
`a022-gpu-tau2-lift-calibration-r022`: n=3/P4096/R64, n=5/P8192/R64, and
n=7/P16384/R64, all at tau=2 and L=16. It includes each complete raw JSONL,
sidecar, strict Decimal-80 aggregate, exact-comparison report, and monitor
receipt. Every raw file ends in `run_end`; every monitor ends with
`returncode: 0`.

The root-audited kernel source is the parent package's `../tau2_lift.swift`
(SHA-256 `e268035742b1c8eb7b0e8c366e3b8e393054352ed0755cd523770e1d8d7d9c3a`),
and the strict aggregator is `../aggregate_tau2_lift.py` (SHA-256
`a0d9e3a602dfcab33735c750104b63e9a5949d0e87cc89f0f8a2c6887b385307`).
For n=3,5,7 the strict estimates of Q are respectively 0.4896687010,
0.3834728035, and 0.3443236389. These are finite simulation diagnostics, not
confidence statements or asymptotic proof.

Relative replay examples:

```sh
python3 -B ../aggregate_tau2_lift.py raw/n3-p4096-r64-l16.jsonl \
  --sidecar receipts/n3-p4096-r64-l16.sidecar.json --source ../tau2_lift.swift
python3 -B compare_exact_relative.py
```

`compare_exact_relative.py` is a clearly named relocation of the retained
historical comparator; it differs only in certificate paths and its numeric
output agrees with `reports/exact-comparison.json` after those paths are
normalized. The corrected report records n=3 Q absolute error
`1.27090000514e-4`; the earlier `1.16e-4` was the normalizer error, not Q.
