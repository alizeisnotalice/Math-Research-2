# n=16 固定根恢复批

原 48 个 replicate 因重启中断后的恢复以原 seed 补齐至 64。`canonical-64-ledger.csv`
汇总完整账本；前 48 行来自 `reference/original-reps0-47.csv`，后 16 行各由 `raw/`
中的单行 CSV 与资源回执重放。`post_recovery_audit.py` 已改为只读此包内相对输入。

审计验证 seed 公式、64 个连续 replicate ID、dead=0、正有限 normalizer、终端 hash、
16 个恢复 raw/资源回执、以及退出状态。归档的 ESS/normalizer 只是有限权重诊断，
不能认证混合、无偏性或原条件律。
