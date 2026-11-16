# r012 来源索引与归档白名单

本包是一个独立、紧凑的原条件律 \(\tau\)-SMC 审计包。下列文件是可归档白名单；
二进制、Swift/Metal 模块缓存、编译日志和临时输出均不在其中。

| 文件 | 角色及证据范围 | SHA-256 |
|---|---|---|
| `DESIGN.md` | 模型、FK 恒等式、一般 offspring-均值无偏归纳及限制 | 见 `HASHES.sha256` |
| `README.md` | 已运行参数、结果、计时与限制 | 同上 |
| `original_tau_smc_audit.py` | CPU 整数穷举与受限树交叉核验 n=2,3 | 同上 |
| `tau_smc_audit_final.json` | 前项的实际保存输出 | 同上 |
| `source/tau_smc_integer_n7_used.swift` | 合法 n=7 输出所用 Swift/Metal 源，hash `f5210b8c…` | 同上 |
| `source/tau_smc_integer_n8_timed_used.swift` | n=8 timed raw 所用源，hash `3de43d78…`；其旧字段名须按 README 更正解释 | 同上 |
| `source/tau_smc_integer_guarded_future.swift` | 未替换已运行生产的后续 guard/术语修订源；仅 smoke 验证 | 同上 |
| `source/summarize_tau.py` | 合并各 \(\tau\) 的零权重保留摘要；旧字段转换为单步直接父诊断 | 同上 |
| `source/normalize_dead_logz.py` | 仅将显式 `dead=true,Z=U=0` 的 `logZ:-inf` 转为 JSON `null`，否则失败 | 同上 |
| `source/history_debug_audit.py` | CPU 完整标签小池转移/峰验证；非 GPU 样本审计 | 同上 |
| `data/n2-valid-smoke.jsonl` | GPU n=2 smoke 原始记录 | 同上 |
| `data/n3-valid-p512-r64.jsonl`, `data/n3-valid-summary.jsonl` | GPU n=3 有限重复原始记录及摘要 | 同上 |
| `data/n7-p16384-r8.jsonl` | 保留的历史 raw；含非 JSON `-inf`，不可直接作 JSON 输入 | 同上 |
| `data/n7-p16384-r8-normalized.jsonl`, `data/n7-summary-direct-parent.jsonl` | 受限转换后的 n=7 合法 JSON 与纠正术语摘要 | 同上 |
| `data/n8-p32768-r16-timed.jsonl`, `data/n8-summary-direct-parent.jsonl` | n=8 timed 原始记录及纠正术语摘要 | 同上 |
| `data/history-debug-n4-p64-r2.json` | CPU 完整标签 debug 的26个保存轨迹 | 同上 |

## 独立小规模锚点

`original_tau_smc_audit.py` 对全部均匀下降标号排列直接计数，并独立以受限前缀
树计数固定 \(\tau\) 的事件；输出记录所有非零层逐项一致：

| n | \(A_2\) | \(Q=\mathbb E[d_E^2/n^2\mid K=2]\) |
|---:|---:|---:|
| 2 | 8 | \(1/2\) |
| 3 | 143808 | \(1100/2247\) |

这是 n=2,3 的有限精确核验，不外推为 n→∞ 结论。GPU n=2 smoke 重现
\(Z_2=1/3,U_2=1/6\)；GPU n=3 的有限重复只作数值比较。

## 仍然存在的限制

1. 理论一阶无偏性以理想随机数、正确 mutation 核和每父 offspring 条件期望为假设；
   实际有限 PRNG 和固定 seed 只给数值证据。
2. 任意有限 \(P,R\) 的 \(\sum\widehat U/\sum\widehat Z\) 不是无偏条件矩估计；
   没有给出 n=8 的严谨置信区间或精确基准。
3. 旧 `min_distinct_ancestors` 及新摘要的 `min_one_step_distinct_parents` 都不是跨代
   祖先数；只能作一步重采样退化诊断。
4. n=8 原始生产没有完整标号，CPU history debug 不能替代实际 GPU 轨迹的独立峰审计。
5. GPU 时间戳、同步提交等待累计和全进程墙钟来自不同聚合层次；其比值不衡量整机 AGX
   利用率，也不可相加。
