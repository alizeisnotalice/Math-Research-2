# 第 28 轮待整合：d=8 分数前缀 GPU 批次

本包保存完成的 `m=36,n=68,d=8` 分数前缀浮点 GPU 探索。低阶模型允许每个 (i<8) 的 pair 为 `(0,0),(1,0),(0,1)`；只保留低卷积合法并通过声明的严格必要过滤条件者。它不是严格排除或精确反例。

独立验证器不导入 worker 或 prefix module，直接枚举并检查所有 2,187 个 ternary 前缀：531 个低端合法，421 个通过必要过滤；对应 SHA-256 为 `2dd97282cba953461ce00dc76114364b56b3b5b33e0b0a4c4b211a727f788c2e` 和 `ac2b52c45594923d316b7709ae0f9610b68fea87c8b22cc4019f4adc0b26a701`，与 CPU audit/config 一致。每 cohort 的 32,768 starts 覆盖全 421 池。

账本为 5 restart 快照 80 行加 final complete cohort 16 行，共 96 个互异保存样本、partial=0。对每行以 Fraction 重卷积并校验 assigned/derived prefix、低卷积、d-cross、端点、非负及 margin。Fraction 与保存 residual 的最大差 `4.4e-17`；最小保存 Fraction residual `0.00211411714553833`，故无精确乘积。

`parent-exit-verification.json` 记录 session 66912 rc0、parent 4557/child 4561 不存在。退出前 allocator 值只记录于 `release.json`，不代表最终常驻内存。

保留 config/checkpoint/CPU prefix audit/release/queue/exit、5 restart、最终 worker/driver/module与被引用旧 worker、独立 verifier/receipt；排除 state、日志和 result。复现：

```sh
python3 scripts/verify_prefix_fraction.py \
  --run-dir run \
  --worker scripts/gpu_factor_search_prefix.py \
  --prefix-module scripts/fractional_prefixes.py \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/prefix-d8-recheck.json
```
