# r022 tau=2 lift implementation and aggregation certificate

Sources: `a022-tau2-lift-root-audit-r022`, `a022-gpu-tau2-lift-r022`, and
`a022-tau2-lift-aggregate-r022`.  This package contains only small kernel
audit receipts and the final strict aggregator
`aggregate_tau2_lift.py` (SHA-256
`a0d9e3a602dfcab33735c750104b63e9a5949d0e87cc89f0f8a2c6887b385307`).
Large n=32 trace JSON and every binary are deliberately excluded.

The retained n=3/P17/R2 raw integration input and sidecars yield the small
finite diagnostic recorded in `n3-p17-r2-strict-summary.json`; it is not a
production run or asymptotic inference.  Reproduce with
`python3 -B aggregate_tau2_lift.py raw/n3-p17-r2-new.jsonl --sidecar
n3-p17-r2-verified.sidecar.json --source tau2_lift.swift`, then
`python3 -B test_tau2_lift_aggregate.py`.  New calibration work remains out
of scope and is not represented here.

The completed r022 n=3/5/7 tau=2 calibration is in `calibration/`; it has
complete raw/sidecar/telemetry/strict-summary receipts and a separate
`calibration/ARCHIVE.md`. It excludes all unfinished n=12 geodesic work.
