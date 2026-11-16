# 第 26 轮：boundary 自由系数 d=1 批次归档

本包是 `m=12,n=24,d=1`、`margin=0.001`、32,768 starts 的单完整 cohort 浮点投影搜索。它不是严格证书，也不是精确反例。

独立脚本把保存 JSON 十进制值转为精确 `Fraction` 后重新卷积 final checkpoint 的 16 个保存样本：`restart_files=0`、`complete_final_cohort_saved_rows=16`、`partial_saved_rows=0`。端点为 1、非负、指定分数位落在实际 `[0.001,0.999]` margin；16 对均含分数系数、精确二元对为 0。保存 float64 residual 与 Fraction 卷积的最大差为 `8.719249979494215e-17`；样本中的最小 Fraction residual 是 `0.0015814881152234545`，故它明显不是精确乘积。

`parent-exit-verification.json` 记录 session 6195 rc0/reaped、PID 74196 不存在。`release.json` 的 `29098752` bytes 是退出前 allocator 指标，进程不存在才是最终释放证据。

最小包包含 config、final checkpoint、release、parent exit receipt、冻结 worker、独立 Fraction verifier 和实际 receipt；明确不含 `state.pt`、`result.json`、日志或任何 round027 产物。没有 queue/driver，因为本批次直接启动 worker。

```sh
python3 scripts/verify_free_factor_fraction.py \
  --run-dir run \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/boundary-free026-recheck.json
```

worker SHA-256 为 `2080fcc0fbfb72040f736606878d02e8f9f307f9f1b18abb43a8514a4ee8a4cf`，与 config 一致。
