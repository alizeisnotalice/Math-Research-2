# 第 26 轮：自由系数 d=2 GPU 批次归档

本包保存已完成的 `m=24,n=64,d=2` 浮点投影搜索。它是指定低阶支撑、固定次数、有限随机 starts 的数值探索；不构成不可行性证明或精确反例。

## 独立 Fraction 审计

独立脚本不导入 GPU 搜索器，将每个保存 JSON 系数的显示十进制值转换为精确 `Fraction` 后卷积。共 5 个轮换快照（80 行）和一个**完整最终 cohort**（16 行），共 96 个互异保存样本；`partial_saved_rows=0`。这 6 个 cohort 由 `12000 / 2000 = 6` 得到，不能把最终 16 行误列为 partial。

所有保存样本端点为 1，两个低阶前缀系数 (1,ldots,d-1) 均为 0，其余保存系数均在 `[0,1]`；指定 (d=2) 系数位于 `[0.05,0.95]`。指定互补和相对 1 的最大 JSON-float32 偏差为 `2.980232244e-08`。96 对都包含分数系数，精确二元因子对为 0。Fraction 卷积的最大 residual 与保存 float64 诊断相差至多 `2.2293649167901823e-16`。

这只验证保存 JSON 及其浮点诊断的一致性，不能验证未保存张量，不能将十进制 Fraction 检查解释为代数精确证书。

## 退出与内存

`parent-exit-verification.json` 记录 session 48011 reaped、exit code 0、parent 30216 和 child 30217 均不存在。`release.json` 记录 cleanup/exit 意图，退出前 allocator 值 `35226368` 不代表进程退出后的常驻内存；已消失的进程才是最终释放证据。

## 最小文件与复现

包含 config、final checkpoint、5 restart、release、queue、parent exit receipt、冻结 worker/driver、独立 Fraction verifier 和实际 receipt。明确排除 `state.pt`、日志和 `result.json`（其内容与 checkpoint 重复）。

```sh
python3 scripts/verify_free_factor_fraction.py \
  --run-dir run \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/free026-d2-recheck.json
```

worker SHA-256 为 `2080fcc0fbfb72040f736606878d02e8f9f307f9f1b18abb43a8514a4ee8a4cf`，与 config/queue 一致。
