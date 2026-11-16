# 第 25 轮：自由系数 GPU 批次归档

本包是 `m=20,n=80,d=1` 的停止后数值账本。搜索器在每 2,000 步保存一次 16 个最优数值样本；轮换到 40,000 步后主线程中断并回收子进程。它是有界、投影的浮点数值探索，不能证明不可行性，也不构造精确反例。

## 实际分数重算

独立脚本不导入 GPU 搜索器，将 JSON 中显示的每个浮点十进制数解释为精确 `Fraction`，并重新卷积。结果：20 个轮换快照给出 320 个保存行，最终 checkpoint 是下一轮换的 16 个部分行；共 336 个互异 JSON 样本。所有保存样本端点为 1、其余系数在 `[0,1]`、指定 (d=1) 系数处于 `[0.05,0.95]`；指定互补和的最大 float32 输出偏差为 `2.980232244e-08`。所有 336 对均含分数系数，乘积精确二元的对为 0。Fraction 卷积与已保存 float64 residual 的最大差为 `3.109767861997025e-16`。

`receipt/m20n80d1-fraction-receipt.json` 是本轮实际重跑收据。这里的“精确”仅指对保存的 JSON 十进制值做有理算术，不能回溯 GPU 未保存张量或把浮点候选升级为代数证书。

## 退出和内存记录

`release.json` 显示 `cleanup_called=true`、`process_will_exit=true`，但其退出前 `remaining_tensor_bytes_before_exit=39879424`。审计时原 parent `62159` 与 child `62160` 均不存在，故进程退出而非该瞬时 allocator 指标才是最终释放证据。`queue.json` 保留原驱动状态 `running` 和 job 的 `driver_interrupted_after_owned_child_recovery`/returncode 1；这与已确认的进程消失不应混同为仍在运行。

## 内容和复现

最小白名单是 `run/config.json`、`run/checkpoint.json`、`run/restarts/*.json`（20 份）、`run/release.json`、`run/queue.json`、`scripts/*.py`、receipt、清单与本 README。明确不含 `state.pt`、驱动日志或任何新运行产物。

```sh
python3 scripts/verify_free_factor_fraction.py \
  --run-dir run \
  --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/free025-recheck.json
```

`config.json` 与 `queue.json` 中记录的 worker SHA-256，和冻结 worker 源文件的实际 SHA-256，均为 `2080fcc0fbfb72040f736606878d02e8f9f307f9f1b18abb43a8514a4ee8a4cf`。
