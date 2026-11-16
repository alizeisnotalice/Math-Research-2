# 第 26 轮：自由系数 d=3 GPU 批次归档

本包保存已完成的 `m=28,n=60,d=3` 浮点投影搜索。它是指定低阶支撑、固定次数和有限随机 starts 的数值探索；不构成不可行性证明或精确反例。

## 独立 Fraction 审计

独立脚本不导入 GPU 搜索器，将保存 JSON 的显示十进制系数解释为精确 `Fraction` 并重新卷积。5 个轮换快照保存 80 行，final checkpoint 是第 6 个完整 cohort 的 16 行；合计 96 个互异样本，`partial_saved_rows=0`。

所有保存样本端点为 1，两个强制低阶前缀系数 (1,2) 均为 0，其他系数在 `[0,1]`，指定 (d=3) 系数位于 `[0.05,0.95]`。指定互补和相对 1 的最大 JSON-float32 偏差为 `2.980232243e-08`。96 对均有分数系数、精确二元因子对为 0。Fraction 卷积 residual 与保存 float64 诊断的最大差为 `2.3745431032782654e-16`。

这里的 Fraction 只审计保存的 JSON 十进制值，无法核验未保存张量，不能升级为代数严格证书或一般排除结论。

## 退出与最小文件

`parent-exit-verification.json` 记录 session 52410 reaped/rc0，parent 45768、child 45769 均已消失。`release.json` 的 `34276096` bytes 是进程退出前 allocator 指标；进程不存在才是最终释放证据。

保留 config、final checkpoint、5 restart、release、queue、parent exit receipt、冻结 worker/driver、Fraction verifier 和 receipt。明确排除 `state.pt`、日志、`result.json`，且未读取正在执行的 d1 boundary 批次。

```sh
python3 scripts/verify_free_factor_fraction.py \
  --run-dir run \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/free026-d3-recheck.json
```

worker SHA-256 为 `2080fcc0fbfb72040f736606878d02e8f9f307f9f1b18abb43a8514a4ee8a4cf`，与 config/queue 一致。
