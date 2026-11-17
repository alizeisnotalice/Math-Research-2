# Replay

Run from the final r027 certificate directory. The proof and production-audit packages refer relatively to r026.

```sh
# Set A022_REPLAY_DIR to an existing temporary directory for generated outputs.
: "${A022_REPLAY_DIR:?Set an existing temporary output directory}"
python3 -B r027-monotone-n7/audit_scan.py
python3 -B r027-monotone-proof/audit_adjacent_turn_order.py
python3 -B r027-monotone-proof/audit_contraction_turn_edges.py
python3 -B r027-monotone-proof/audit_four_word_turn_convexity.py
python3 -B r027-merge-rb-production-audit/test_production_audit.py
python3 -B r027-l64-recovery/cone/audit_cone.py \
  r027-l64-recovery/canonical/n24-move64-canonical-eps-1-2.json
(
  cd r027-global-normalizer
  python3 -B punctured_normalizer_table.py 3 --output "$A022_REPLAY_DIR/n3.json"
  python3 -B verify_rank_one_deletion_n3.py
)
(
  cd r027-completion-layer
  python3 -B pratt_verify.py
  python3 -B prime_audit.py
  python3 -B crt_verify.py modular-runs
  python3 -B n8-remote/audit_receipts_portable.py
  python3 -B n8-independent/verify_n8.py pratt-certificates.json \
    n8-remote/remote-results/n8-p0.json n8-remote/remote-results/n8-p1.json \
    n8-remote/remote-results/n8-p2.json n8-remote/remote-results/n8-p3.json
)
python3 -B r027-n12-merge-rb/source/validate_n12_trace.py \
  r027-n12-merge-rb/trace/n12-p17-r1-l16.jsonl \
  r027-n12-merge-rb/trace/n12-p17-r1-l16.jsonl.trace.json \
  ../r026/r026-merge-rb-calibration/source/tau2_merge_rb_calibration_r026.swift
python3 -B r027-n12-merge-rb/source/aggregate_production_r026.py \
  r027-n12-merge-rb/raw/n12-p131072-r32-l16.jsonl \
  ../r026/r026-merge-rb-calibration/source/tau2_merge_rb_calibration_r026.swift
```

The completion-layer source can be compiled for an independent finite pass;
its n=8 mode requires an explicit flag and is not a replay target for this
round. The L64 cone check is an exact rational necessary-condition check of
the stored plug-in summaries, not a confidence interval or an accuracy proof.
