# 第47轮三次 ANF 五维平坦的双机整数枚举

2026-11-18，本机是 12 逻辑核的 Mac mini，远端是 15 逻辑核的 MacBook Pro。两机均用 arm64 macOS 和 Apple Clang 21.0.0。对每个九变量三次 ANF，枚举全部 3,309,747 个五维线性子空间及各自 16 个平移，共 52,955,952 个五维仿射平坦。

两个完整枚举独立排列 RREF 主元型与非主元参数。第一种算法直接求 ANF 真值，按 64 样本将真值打包为位平面，以按位交、并判断单色。第二种以 Möbius 变换生成逐样本的 512 点位集，以 popcount 完整统计每个限制的 Hamming 重量。两种算法的单色计数与重量端点逐样本一致。

首批 64 个样本采用种子 202611180958，其中 25 个在任何五维仿射平坦上都不恒定。Python 独立验证高斯乘积和递推、分片区间、5,280 个试跑平坦、77 个分片正见证，以及每个重量直方图的前四阶阶乘矩。输入一行对应一个样本，130 个系数按单项式掩码升序排列。精确计数、见证和直方图见 combined_result.json。

扩展批以种子 202611181006 产生 32 组共 2,048 个新样本。每种算法独立处理 108,453,789,696 个样本–平坦对，发现 849 个样本没有单色五维平坦。验证全部扩展样本的前三阶关联矩与两个算法端点一致性，并逐项复核 768 份跨机器重复分片。扩展样本逐个计数、输入哈希和负载摘要见 extended_summary.json。

## 负载和退出

两个基准枚举均由两机同步分片，MacBook 运行 0–14 号分片，Mac mini 运行 15–26 号分片。首批直接位平面程序两机墙钟分别为 0.0496 和 0.0358 秒。相同样本的独立全重量谱两机墙钟分别为 0.5928 和 0.5483 秒。

扩展批最初的 MacBook SSH 连接在 banner 阶段超时，当时 Mac mini 完成了其互补分片。严格主机校验的 SSH 恢复后，MacBook 完成剩余 15 片，Mac mini 同步跨机器复算其中 12 片，二者重叠 19.841 秒。两机系统稳态采样中位数均为 100%。按进程累计 CPU 时间及墙钟计算，MacBook 和同步复算的 Mac mini 在重量谱计算阶段平均分别使用 6.98 与 5.12 个等效核，采样窗口峰值分别为 11.63 和 10.62 个等效核。Mac mini 首轮互补分片的采样峰值为 11.78 核。两机上另有其他研究计算，均未干预。

所有枚举进程正常退出。逐机检查无本批进程残留，扩展远端 131 个输出文件的 SHA-256 全通过，三次扩展运行的标准错误均为空。

## 复现

在各自工作目录生成二进制和输入。保留断言。本机的共享编译缓存如下，远端使用其已有的 "$HOME/CodexRemoteMath/cache/clang"。

    r047_module_cache=/Volumes/data/CodexTemp/cache/shared-clang-modules
    clang -O3 -std=c17 -pthread -fmodules-cache-path="$r047_module_cache" flat_census.c -o r047_flat_census
    clang -O3 -std=c17 -pthread -fmodules-cache-path="$r047_module_cache" weight_census.c -o r047_weight_census
    python3 generate_samples.py anf_samples.txt
    python3 generate_extended.py extended_input 32

以下命令中的 0 为立即运行。同步实测时将两机对应的 0 换为同一未来 UTC Unix 时间戳。

    ./r047_flat_census anf_samples.txt 10000 0 1 12 0 pilot.json pilot_usage.csv
    ./r047_flat_census anf_samples.txt 27 15 12 12 0 mini_shards.json mini_usage.csv
    ./r047_flat_census anf_samples.txt 27 0 15 15 0 macbook_shards.json macbook_usage.csv
    ./r047_weight_census anf_samples.txt 27 15 12 12 0 mini_weights.json mini_weights_usage.csv
    ./r047_weight_census anf_samples.txt 27 0 15 15 0 macbook_weights.json macbook_weights_usage.csv
    bash run_extended.sh mini 15 12 0
    bash run_extended.sh macbook 0 15 0
    bash run_extended.sh mini_cross 0 12 0

将两机的输出合到同一目录后进行独立核验。

    python3 verify_and_aggregate.py . combined_result.json
    python3 verify_extended.py . extended_summary.json

计时和负载会随运行变化。证书涵盖指定样本和全部五维平坦。
