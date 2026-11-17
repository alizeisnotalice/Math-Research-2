# 实际重放记录

归档前已经实际通过的可携带重放如下（均在各自目录中执行）：

```sh
python3 -B r025-fixed-root-n11/verify_archive.py
python3 -B r025-strip-entropy/audit_strip_entropy_n3.py
python3 -B r025-rank-bias/verify_rankbiased_prefix.py
python3 -B r025-grid-reflection/audit_n3_root_moments.py
python3 -B r025-grid-geometry/audit_toprow_singlepeak_probability.py
python3 -B r025-grid-geometry/audit_two_row_completion_n5.py
python3 -B r025-genealogy/diagnose_n24.py
python3 -B r025-normalizer-tail/audit_n3_normalizer_tail.py
python3 -B r025-normalizer-tail/audit_n24_finite_diagnostic.py
python3 -B r025-normalizer-tail/certify_plugin_error_floor.py
python3 -B r025-rank-bias/verify_rankbiased_trace.py \
  r025-rank-bias/n3-p17-r1-l4-trace.jsonl \
  r025-rank-bias/n3-p17-r1-l4-trace.jsonl.trace.json
gunzip -t r025-macbook-gpu/n16-distance-p262144-r128-l16.jsonl.gz
python3 -B r025-local-n24/aggregate_tau2_distance.py \
  r025-local-n24/n24-p131072-r64-l16.jsonl \
  --sidecar r025-local-n24/n24-p131072-r64-l16.sidecar.json \
  --source r025-local-n24/tau2_distance_terminal_r023.swift --epsilon 1/4
```

`diagnose_n24.py` intentionally reads the two named completed raw inputs from
the same certificate package (`r025-local-n24/` and `r025-macbook-gpu/`) so
that it can compare the two methods; its archived source and replay stdout are
retained, and the report states that its conclusion is a finite diagnostic.

The rank-biased GPU trace was independently replayed before archiving; its
small raw, GPU-side trace payload, and `n3-p17-r1-l4-trace.verify.json` are
all retained.  No run above proves mixing, a confidence interval, or an
asymptotic probability statement.
