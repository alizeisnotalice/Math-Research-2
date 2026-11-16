# r022 GPU terminal-geodesic certificate

Sources: `a022-gpu-geodesic-prefix-r022`, `a022-geodesic-root-audit-r022`,
and the completed portion of `a022-geodesic-calibration-r022`.  The final
source is `geodesic_prefix_r022.swift` (SHA-256 `97b0cb005a85d7ec1dc04800bf1bf0cd338ff1046e096bb4991c1ff7d763a766`).
It is accompanied by small independent-oracle receipts, not large traces or
binaries.

`calibration/` contains only natural-ended n=3/P4096/R64, n=5/P8192/R64,
and n=7/P16384/R64 runs, each with raw JSONL, sidecar, and Decimal summary.
The currently running n=12/P131072/R128 calibration and any partial output
are explicitly excluded.  This package contains finite implementation and
calibration evidence only.  Relative replay commands: `python3 -B
oracle/n3_onepeak_oracle.py`; `python3 -B summarize_onepeak_geodesic.py
calibration/n3-p4096-r64-l16.jsonl`; and `python3 -B
oracle/verify_geodesic_trace.py TRACE.json RAW.jsonl` when a separately
generated trace is supplied.  No GPU job is started by this archive.
