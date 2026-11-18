# 第45轮复现

PKG 设为本目录。OUT 设为挂载数据盘上的新目录。Python 3.15 和 Apple clang 21 已在本机运行。

    mkdir "$OUT"
    env CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang++ -std=c++20 -O3 -Wall -Wextra -Werror "$PKG/rootwise/cut_region_rootwise.cpp" -o "$OUT/rootwise"
    export A022_R045_REGION_BINARY="$OUT/rootwise"
    python3 -B "$PKG/audit_rounding_capacity.py"
    python3 -B "$PKG/rootwise/verify_vectors.py"
    python3 -B "$PKG/rootwise/audit_rootwise.py" "$OUT/n5-gate.json"
    python3 -B "$PKG/rootwise/run_rootwise_n7.py" "$OUT/n7" --workers 12

上面四项核验均在最终公开目录执行。n7 全扫描为 644 个物理计算，使用反射恢复 1220 条曲线。n8 MacBook 的整批为 78 项，实际采用3个 worker，单进程峰值约837 MB。公开代码只改输入路径和执行文件参数。n8 可按下列命令重放，旧汇总用于校验其16条子族。

    python3 -B "$PKG/rootwise/run_parity_n8.py" "$OUT/n8" --workers 3 --old-summary "$PKG/../r044/n8-summary.json"

隔行两点计数分别编译两个源文件，定向容斥为主，直接子集递推为小规模独立交叉核验。定向容斥峰值状态受驱动预算及程序硬上限双重约束。

    env CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang++ -std=c++20 -O3 -Wall -Wextra -Werror "$PKG/domino/count.cpp" -o "$OUT/count"
    env CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pedantic "$PKG/domino/count_ie.cpp" -o "$OUT/count_ie"
    python3 -B "$PKG/domino/run.py" "$OUT/domino.json" --workers 8 --counter "$OUT/count_ie" --direct "$OUT/count"

原始计数及资源收据在各机工作目录保留。公开汇总包含稀疏根向量供独立重聚合，无需读取随机样本。
