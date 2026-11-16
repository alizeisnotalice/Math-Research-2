# 第 27 轮：d=7 分数前缀 GPU 批次归档

本包保存完成的 `m=32,n=72,d=7` 分数前缀模型。它不同于旧的“全零低前缀”模型：对每个 (i<7)，((a_i,b_i)) 可为 `(0,0),(1,0),(0,1)`，随后只保留低卷积合法且通过已声明严格必要过滤的 prefix。GPU 仍是投影浮点数值探索，不构成严格排除或精确反例。

## 独立精确 prefix/Fraction 审计

独立验证器不导入 worker 或 prefix 模块，重建整个有限池：原始 ternary 代码 729，合法低 prefix 227，严格必要条件幸存 prefix 187。合法/幸存列表 SHA-256 分别为 `a5654cf6dd4321be38e54055db807f5a663c51ee310391e398d608c1fcd3f5bf` 与 `d12e7e29e75dd798c573307843a80c2d620b3f9787a27d907506a567b44c252d`，与冻结 CPU audit/config 一致；每 cohort 的 16,384 starts 覆盖全 187 池。

保存账本是 5 restart 快照 80 行加 complete final checkpoint 16 行，合 96 个互异样本、`partial_rows=0`。对每个保存样本，验证器从实际系数精确重建 assigned/derived prefix、低卷积二元性、(d) 交叉项为零、端点/非负/margin，并以 `Fraction` 重卷积。保存 residual 与 Fraction residual 的最大差为 `9.775971246752306e-19`；最小保存 Fraction residual 为 `0.002617187099531293`，故不是精确乘积。

“187 个幸存”只代表此 worker 的两个必要过滤条件均未触发；它不是一般可延拓性或完整非负因子结论。Fraction 算术也只作用于保存 JSON 十进制值，不能回溯未保存 GPU 张量。

## 退出、最小白名单与复现

`parent-exit-verification.json` 记录 session 55630 rc0/reaped、parent 79529 与 child 79530 均不存在。`release.json` 的 `39584512` bytes 是退出前 allocator 值，不能当作退出后内存占用。

包保留 config、checkpoint、5 restart、CPU prefix audit、release、queue、exit receipt、最终 worker/driver/prefix module、被引用旧 worker、独立 verifier 与 receipt。明确不含 `state.pt`、worker log 或 `result.json`。

```sh
python3 scripts/verify_prefix_fraction.py \
  --run-dir run \
  --worker scripts/gpu_factor_search_prefix.py \
  --prefix-module scripts/fractional_prefixes.py \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/prefix-d7-recheck.json
```
