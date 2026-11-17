# r025 rank-biased ordered-prefix move

`tau2_rankbiased_prefix_r025.swift` 从冻结 source
`01b77bca1735f07293eae5ffeea0326810fdc5ea785c7a9e99d9030cba9b0702`
复制后只修改 Metal `move` 核和输出 metadata/trace implementation 字段。
新 metadata 明确标记 rank 选择规则、旧 source 哈希和适用范围。冻结 r023
聚合器显式要求 `implementation=r023-tau2-terminal-L-distance-histogram` 和
`base_source_sha256`，因此会拒绝此输出。未来若汇总 r025 输出，必须新建并审计
一个同时验证 r025 implementation、`frozen_source_sha256` 和 rank 规则的入口；
不能放宽旧聚合器的 guard。

先运行精确 CPU 核验：

```sh
python3 verify_rankbiased_prefix.py
```

本目录尚未运行 Metal/GPU，也没有 SMC 生产结果。未来的小 GPU trace 必须逐proposal
对照此规则和保存的 fixed-seed CPU reference，随后才能作为实现审计；即使通过也只
证明这一前缀 move 的实现，不证明混合或完整 SMC 结论。

2026-11-17 的实际本机核验以 `Fraction` 完成：\(P_2\) 的所有 rank 1--4 和
\(P_3\) 的所有 rank 1--5 前缀，分别对本规则及一组任意固定非均匀 rank 权重
逐行核验行和、支持和详细平衡。输出
`exact-prefix-balance.json` SHA-256 为
`51eaf1b22c7470cf8c9b43991309266d762c4b44a8f6ac53581591ebf291b321`。
Swift 宿主程序仅编译、未运行 Metal：source SHA-256
`60746225d734d03a5567c0b2e45581e89cd8a39599605f6e1a4397b612e23334`，
本机二进制 SHA-256
`94cd9a14484b29d08662bbf7a22193b789c20503e9411cdd0cea7c4b09bb779d`。

随后完成一个独立的最小 GPU 端到端审计：\(n=3\)、\(P=17\)、一个 \(\tau=2\)
重复、`strat`、4 个 move、domain `9e25090000000001`，并启用完整 trace。Swift
运行时 `makeLibrary` 实际编译了 Metal 字符串，进程以退出码 0 结束，墙钟 3.77 秒，
峰值 footprint 40,239,824 bytes、零 swap/零 signal。独立 Python 使用相同的
SplitMix64、fully-adapted resampling、rank-biased proposal 与前缀状态重放 9 个
rank 的全部 17 条轨迹，逐字段匹配 `allowed`、CDF、parent、前缀和 bitset state。
它还实际遇到并检查了根–非根活跃交换、根–外部替换及拒绝：前者 146 次提议
（48 接受、98 拒绝），后者 217 次（134 接受、83 拒绝）。这验证有限实现路径，
仍不验证混合或 SMC 精度。

trace 的 9 个 rank 行不是 9 次完整前缀 move：在 \(n=3\) 时，move 只发生在
rank 1 至 8 的 resample/append 后，rank 9 只执行 terminal evaluator。因此该
试验实际重放 \(17\times8\times4=544\) 次 move 提议；第 9 行仅用于完整的
\(N-1\) 前缀和终端整数核验。

源码中 lift/distance 字符串附近的“frozen SMC, RNG and MH ... untouched”是从
冻结 r023 终端模块继承的历史说明，描述该模块的来源版本；它不描述当前 r025
完整程序。当前 r025 的主 Metal `move` 核确已按本目录所述修改，且该修改正是
本 CPU/GPU trace 审计的对象。
