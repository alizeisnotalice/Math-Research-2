# 第14轮：n8/n9 原条件律矩与一峰归一化

这里归档完整、已经结束的有限批次。`root-audit.json` 由主 Agent 的
`../../computation/audit_r014_gpu.py` 以 Decimal 60 位重新核算；整数种子、完整出生秩、零质量重复、run_end、Metal 字节一致性均检查。
所有标准误为跨重复经验诊断，不是置信认证，也没有把粒子当独立样本。

## 本轮新双峰批次

每个出生秩 tau=2,...,n² 都有16个重复，整数分层重采样。

| n | 粒子数 | seed域（十六进制） | Q | 经验delta SE | 墙钟秒 |
|---:|---:|---|---:|---:|---:|
|8|65536|6100000000000000|.3345058002|.00076231|112.648|
|8|131072|6200000000000000|.3349364452|.00051492|129.806|
|9|131072|6400000000000000|.3288934808|.00100218|211.390|
|9|131072|6500000000000000|.3301308719|.00084399|264.750|

实际种子为 `(0xA022000000000000 + 域 + n*1000000 + tau*10000 + rep) mod 2^64`。
每文件都保留势质量 Z、空间矩 U、死端、GPU时间戳及 host submit+wait 时间。
大批次生产不保存完整标号，不声称对所有生产粒子独立判峰。

`*-monitor.csv` 为每秒整机非独占 AGX 与目标进程 RSS。n8/P65536 的 AGX 字段因早期正则失败为空，记为不可用。
其余三批有效读数分别117、191、239个，均值99.658%、99.743%、99.063%；最后一批含退出边缘零读数。
进程 RSS 不代表全部驱动/Metal 统一内存。GPU时间与墙钟不相加，离散整机读数不证明本作业独占连续满载。

## 一峰与同粒子数合并

`onepeak/source/onepeak_smc.swift` 取 tau=n²+1，仅输出一峰 Z。
一峰种子为 `(域 + n*1000000 + rep) mod 2^64`，没有双峰的 A022 前缀或 tau 项。
小格 n2/n3 与 n7 精确 A1 校准；n8/n9 各用 P131072、每域128重复、8100...和8200...两域，均无死端。
身份、未归一化无偏性的理想随机假设、比值限制见 `../../proofs/one-peak-smc-normalizer.md`。

主汇总只合并同一 P=131072 的批次：n8 复用第13轮7600...域32重复，再加本轮6200...域16重复；n9合并本轮两域。
该旧输入通过路径明确列出，不算本轮新运行；n8/P65536不参与合并。

| n | 双峰每tau重复 | 一峰重复 | Q合并 | R=A2/A1合并 | R经验delta SE |
|---:|---:|---:|---:|---:|---:|
|8|48|256|.3349933536|3933.7624|12.0172|
|9|32|256|.3295104077|15568.1813|87.2255|

它们不判定 R 或 Q 的渐近。审计共检查7304个本轮及明确复用记录的互异种子；互异 seed 本身不证明随机独立。

## 复现及冻结范围

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-modules \
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules \
  production_r013.swift -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/tau_smc_r014
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/tau_smc_r014 \
  9 131072 16 2 81 /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/n9-r014.jsonl \
  full strat 6400000000000000
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-modules \
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules \
  onepeak/source/onepeak_smc.swift -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/onepeak_r014
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/onepeak_r014 \
  9 131072 128 8100000000000000 /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/onepeak-n9-r014.jsonl
PYTHONDONTWRITEBYTECODE=1 python3 ../../computation/audit_r014_gpu.py
```

源码、实际版本及哈希在本目录和课题 `computation/manifest-r014.json`。
完整新大掩码核审计另见 `kernel-audit/`：n9/n10/n12、P17，不扩大为生产全粒子审计。
没有归档二进制、缓存、早期中断的6100同seed片段、误用多项式的6300批或重复归档第13轮n7批。
正在进行的更大n生产不属于本轮冻结结果；本轮已完成进程均已退出。
