# 线程局部修正版 CPU SMC

最终有效文件为 b1-L0.csv、b1-L4.csv、b1-L16-rerun.csv；每份 n8,P4096,R32,8线程，完整出生秩2..64，共2016记录。基础 seed 为 0xA0221E00；每记录实际 seed 在 CSV 中。L16 首次输出不完整，未归档。CSV 采用双精度计数权重，非精确计数。

主审发现并排除了两代旧结果：(1) 取模与浮点 CDF；(2) 随后的共享 LAST/REJ 竞态。最终源用 rejection-below、整数 CDF 和 thread_local 临时变量。小例 TSan 未复现旧竞态，不能用于替旧结果辩护。

```sh
clang++ -O3 -std=c++17 -pthread rmsmc.cpp -o /Volumes/data/CodexTemp/T2/a022-r018-b1
/Volumes/data/CodexTemp/T2/a022-r018-b1 8 4096 32 16 8 /Volumes/data/CodexTemp/T2/a022-r018-b1.csv 0xA0221E00
python3 audit_cpu_csv.py b1-L0.csv b1-L4.csv b1-L16-rerun.csv
python3 verify_traces.py b1 3 b1-n3-trace.mh.csv
```

完整轨迹复现：从临时目录设置 RMSMC_TRACE=1，运行 n3 P17 R1 L4 threads1；生成 rmsmc-mh.csv。每个移动的 before/after、接受决定由独立邻接集合判定。主摘要 cpu-threadsafe-audit.json 用 Decimal80 直接读原 CSV，取各秩均值之和的比值。*.final.json 内 mean_rep_ratio 和其 SE 是另一个统计量，不能把该 SE 配给主比值。所有 SE 仅描述性；不认证置信覆盖/混合。

本轮只研究自由边界方格、原图距离、均匀标号条件 K=2。有限数据不证明渐近结论。

复现环境为 Python 3.15、Apple Clang 21、Swift 6.4.2 / Metal。所有编译与输出请放 `/Volumes/data/CodexTemp/T2`，Clang/Swift 模块缓存使用 `/Volumes/data/CodexTemp/cache`。

归档版额外加入参数及 TRACE 单线程入口检查，未改变核。原运行源和归档源哈希见 source-provenance.json。主审实际确认 TRACE 两线程被拒绝；单线程 seed0xA0222000 的 CSV 和4284事件轨迹与原源逐字节一致，见 guard-and-trace-check.json。
