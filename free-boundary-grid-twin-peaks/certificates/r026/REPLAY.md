# 可携带重放

在最终证书目录下执行；`r026-consistency-scale` 与 `r026-merge-rb` 通过同级既有
r024/r025 证书引用其已归档输入，避免重复 raw。若在仓外重放，设置
`A022_CERTIFICATES_ROOT` 为含 r024、r025、r026 的证书根。

```sh
python3 -B r026-row-release-formula/audit_row_release.py
python3 -B r026-row-release-native/verify_runs.py
python3 -B r026-row-release-macbook/verify_n6.py
python3 -B r026-row-release-normalization/exact_surface_table.py
python3 -B r026-row-release-normalization/n7_addendum.py
python3 -B r026-row-release-independent-ie/verify_package.py
python3 -B r026-row-release-independent-ie-n7/audit_results.py
python3 -B r026-monotone-prefix-scan/audit_monotone_prefix.py
python3 -B r026-singleton-slice/singleton_slice_identity.py
python3 -B r026-fixed-root-n16-recovery/post_recovery_audit.py
python3 -B r026-consistency-scale/audit_scale.py
python3 -B r026-merge-rb/audit_merge_rao_blackwell.py
python3 -B r026-merge-rb/verify_merge_rb_trace.py \
  r026-merge-rb/r026-gpu-trace/run-n3-p17-v2.jsonl \
  r026-merge-rb/r026-gpu-trace/run-n3-p17-v2.jsonl.trace.json
python3 -B r026-merge-rb-calibration/audit_derived_trace.py \
  r026-merge-rb-calibration/raw/smoke-n3-p17.jsonl \
  r026-merge-rb-calibration/raw/smoke-n3-p17.jsonl.trace.json \
  r026-merge-rb-calibration/source/tau2_merge_rb_calibration_r026.swift
python3 -B r026-merge-rb-calibration/aggregate_calibration.py \
  r026-merge-rb-calibration/raw/n3-p4096-r64.jsonl \
  r026-merge-rb-calibration/source/tau2_merge_rb_calibration_r026.swift
python3 -B r026-merge-rb-calibration/aggregate_calibration.py \
  r026-merge-rb-calibration/raw/n7-p16384-r64.jsonl \
  r026-merge-rb-calibration/source/tau2_merge_rb_calibration_r026.swift
python3 -B r026-merge-rb-calibration/test_aggregate_rejections.py
python3 -B r026-macbook-fixed-root/audit_remote_results.py
```

这些重放分别认证解析恒等式、有限整数结果、账本或实现路径。它们均不认证 SMC 混合、
置信区间、渐近概率或原题。
