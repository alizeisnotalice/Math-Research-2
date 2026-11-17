# 首次合并生产批聚合审计

本脚本只接受哈希为 `3a9f3c0a8c10409346c25f16f2d82463f928bfaacf7307b4a88a770ec8f69ec6` 的派生 Swift 源，适用 n=3..32 的 tau=2 strat/multi 完整输出。它不接受旧 n24 distance-terminal 源的原始记录。

```sh
python3 -B aggregate_production_r026.py RAW SOURCE
python3 -B test_production_audit.py
```

主线复核发现初版使用字典推导会静默覆盖重复 replicate ID；最终版先校验行数、类型和唯一性，再按 ID 对齐，检查 metadata、种子、normalizer、整数群组和 RB 几何。它以 logZ 计算权重，允许实际浮点指数下溢的 Z=0，保留 log-normalizer。

最终实际核验：n3/n7 同一 raw 的 Q 与已归档校准结果一致；9 种异常输入被拒绝；打乱记录顺序仍得到相同结果。人工 logZ=-1000 下溢案例只测试聚合逻辑，不是 n32 数值实验。每个测试的输入哈希及审计器哈希见 production-audit-tests.json。

输入审计不认证 SMC 的精度、置信区间、混合、方差改善或原题的渐近概率。可选 --exact-dir 仅在存在匹配尺寸的输入时提供比较，否则明确缺失；不使用其他尺寸的默认值。
