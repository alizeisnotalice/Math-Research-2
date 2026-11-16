# 百万粒子 n16 已完成单层批次

两份 JSONL 各128重复：一峰 tau257、domain8700000000000000；双峰 tau2、domain7400000000000000。P1048576，integer stratified parents。源使用 ../r017-bitset/bitset_smc_checked.swift。两个批次已结束；正在运行的完整tau2..256批次不在本轮证书中。

```sh
python3 summarize_logz_deadzero.py n16-onepeak-p1048576-r128-d8700.jsonl n16-doublepeak-tau2-p1048576-r128-d7400.jsonl
```

root-summary.jsonl 为 Decimal log域重聚合，明确区分 alive underflow 与 dead zero。ESS 和最大权重仅为重复权重集中诊断。单个 tau2 层不能作为完整双峰条件律，也不能与一峰平均值直接给全 A2/A1。

本轮只研究自由边界方格、原图距离、均匀标号条件 K=2。有限数据不证明渐近结论。

复现环境为 Python 3.15、Apple Clang 21、Swift 6.4.2 / Metal。所有编译与输出请放 `/Volumes/data/CodexTemp/T2`，Clang/Swift 模块缓存使用 `/Volumes/data/CodexTemp/cache`。
