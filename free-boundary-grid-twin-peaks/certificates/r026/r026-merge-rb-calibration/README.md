# r026 first-merge RB finite calibration

Replay:

    python3 audit_derived_trace.py raw/smoke-n3-p17.jsonl raw/smoke-n3-p17.jsonl.trace.json source/tau2_merge_rb_calibration_r026.swift
    python3 aggregate_calibration.py raw/n3-p4096-r64.jsonl source/tau2_merge_rb_calibration_r026.swift > results/n3-summary.json
    python3 aggregate_calibration.py raw/n7-p16384-r64.jsonl source/tau2_merge_rb_calibration_r026.swift > results/n7-summary.json

The source guard distinguishes the frozen r026 trace source from the observed
calibration derivative. The bundled `exact-top2-d4/` directory is an explicit
baseline parameter; aggregate audit success validates files and exact baseline
matching only, never statistical accuracy.
