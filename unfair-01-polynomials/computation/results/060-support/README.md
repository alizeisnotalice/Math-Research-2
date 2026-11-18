# 第 60 轮严格有限支撑范围

| 每侧最大支撑 | 总次 | 原始支撑对 | 规范代表 | 初筛未判 |
| ---: | ---: | ---: | ---: | ---: |
| 6 | 66 | 23,121,321,333 | 11,115,008,263 | 0 |
| 6 | 67 | 25,624,394,114 | 12,812,391,834 | 0 |
| 9 | 41 | 82,208,894,672 | 41,104,672,088 | 142 |
| 9 | 42 | 131,667,572,169 | 61,073,063,597 | 158 |

六点使用[第59轮的128位并集实现](../059-support/wide-engines/)及本目录的固定数组逐和计数实现。MacBook 九点层使用[第56轮并集实现](../056-support-compact/)及同一固定数组实现。两次运行的分片覆盖和分类逐条相同。六点计 4,115 片，运行在本机各 10 进程，耗时 224.67 和 375.80 秒。九点计 18,357 片，运行在 MacBook 各 13 进程，耗时 110.66 和 1099.95 秒。完成的两批主管和工作进程均已退出。

固定数组实现在 7 组对照中对 19,937,223 个原始对的记录与既有两实现匹配，最高总次数 69。冻结源及摘要位于 compact-independent。以 Clang C++20 -O3 编译同目录里的 support_bitset.cpp 为 support_bitset，再用 run_bitset.py 启动。枚举边界及支撑空位传播规则见[第56轮说明](../056-support-compact/README.md)。

九点层 300 个筛后候选中，272 个直接对角证书、27 个有理正恒等式和[最后一例的单独恒等式](nine41-42/case262.md)全部通过逐项精确复核。默认低阶模板在最后一例未得到证书，随后用第 4、25 行消元得到简短的手工证明。

在本目录调用下列命令可从归档分片重建覆盖报告，并逐条验证九点残余证书。使用 Python 3.13.18 和 SymPy 1.14.0。

    python3 -B ../056-support-compact/compare_complete.py six66-67/production six66-67/independent --totals 66 67 --cap 6 --out /Volumes/data/CodexTemp/T2/a30-six67-check.json
    python3 -B ../056-support-compact/compare_complete.py nine41-42/production nine41-42/independent --totals 41 42 --cap 9 --out /Volumes/data/CodexTemp/T2/a30-nine42-check.json
    python3 -B nine41-42/verify_residual.py
