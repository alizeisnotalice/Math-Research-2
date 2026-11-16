# 第16轮：已完成的n10/n12 GPU批次

归档范围：`production/`仅含P262144、R16、四个完整双峰批次，`onepeak/`含两域各256重复的一峰归一化及n7小校准。
P524288、domain7000新批次仍运行，未混入本轮结果。

内核为既有r013/r014已逐秩核验的五字占据掩码Metal核，支持n≤12。本轮没有使用待审计位集新核。
标准字符串提取（保留首尾换行）的SHA256为
`7474e22aee24049eb8692db789f32dff7b46e5c1c9d31f569622a7278c47da39`；
完整生产Swift源为`cdcd733e41bae122666b03c675123963b7a5ae9a26643dfd458493d156063523`。
主Agent核查双峰、一峰、已有冻结核的Metal字节相同。对计数及几何语义的既有小粒子审计参见r014-gpu，未声称本轮保存了每个生产粒子的完整标号。

## 输入和复现

Apple M6、12 CPU核/12 GPU核、16 GiB；Swift 6.4.2 / Metal，Python 3.15。
双峰n10使用域6600000000000000、6700000000000000，n12用6800000000000000、6900000000000000（均为十六进制参数）。
每域全tau=2,…,n²，每tau16重复，P262144；整数分层选父，所有死端零权重保留。
seed=(0xA022000000000000+domain+n×1000000+tau×10000+rep) mod 2^64。
一峰tau=n²+1，域8300…/8400…，P262144、R256；seed=(domain+n×1000000+rep) mod 2^64。

```sh
swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules -framework Metal production_r013.swift -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r016-production
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r016-production 12 262144 16 2 144 /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r016-n12-repeat.jsonl full strat 6900000000000000
# onepeak目录内的独立host命令见onepeak/REPORT.md。
# 从课题根目录独立核算所有保存原始重复：
python3 computation/audit_r016_gpu.py
```

## 核查和资源

四批均有完整run_end；每批最后两个出生秩是dead，其他重复无非dead Z下溢。最小正Z仍约10⁻⁵¹。
主Agent以Decimal60重新核算8784条含n7校准的记录，实际seed全部互异。
独立worker的按tau统计另存worker-aggregate-262144.json；其中“独立”措辞已修正为不同seed重复、未认证随机独立。
n10的tau73两域差在Z/U约3.67/3.93个经验SE，n12最大约2.66/2.23；这是多层有限诊断，不解释成系统效应或置信检验。

1Hz monitor的整机非独占AGX采样（不是专属内核计量）：

| 域 | 有效样本数 | 均值% | 100%样本数 | 墙钟秒 |
|---|---:|---:|---:|---:|
|n10/6600|827|99.6143|813|887.048|
|n10/6700|484|99.7645|469|518.032|
|n12/6800|1172|99.5563|1149|1252.087|
|n12/6900|1171|99.6815|1105|1239.895|

各重复处于autoreleasepool内。run_end后worker核查旧生产进程退出；主Agent再次按确切可执行路径和参数核查本轮生产及一峰进程均已不存在，见process-exit-check.json。
另一个P524288批次在继续，不能将当前总GPU驱动内存归因于已退出批次。进程RSS也不代表全部Metal驱动分配。

全部Q、R和SE只是有限数值证据。理想随机数下FK未归一化无偏，不保证有限PRNG、有限粒子比值或有限重复均值无偏；互异seed也不证明独立。没有严格误差条或渐近结论。
