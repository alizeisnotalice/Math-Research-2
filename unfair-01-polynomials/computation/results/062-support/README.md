# 六点支撑的第 68 至 73 层

两套精确实现按 128 位并集和固定数组逐和计数分别独立重算，每次用本机 10 个进程。

| 总次数 | 原始支撑对 | 约化代表 | 未判 |
| ---: | ---: | ---: | ---: |
| 68 | 30,532,616,148 | 14,694,501,084 | 0 |
| 69 | 33,732,881,448 | 16,866,136,001 | 0 |
| 70 | 39,987,667,156 | 19,265,031,983 | 0 |
| 71 | 44,049,685,135 | 22,025,101,011 | 0 |
| 72 | 51,964,188,379 | 25,058,919,604 | 0 |
| 73 | 57,084,789,769 | 28,542,690,884 | 0 |

第 68、69 层的 5,418 片与第 70、71 层的 7,079 片逐项相同。第 72、73 层的 9,182 片也全部匹配。严格结论为每侧至多 6 点、总次数至多 73 时只能公平。九点的完整已核验界仍是 42。三批的主管和各工作进程均已退出。

    python3 -B ../056-support-compact/compare_complete.py six68-69/production six68-69/independent --totals 68 69 --cap 6 --out /Volumes/data/CodexTemp/T2/six69-check.json
    python3 -B ../056-support-compact/compare_complete.py six70-71/production six70-71/independent --totals 70 71 --cap 6 --out /Volumes/data/CodexTemp/T2/six71-check.json

    python3 -B ../056-support-compact/compare_complete.py six72-73/production six72-73/independent --totals 72 73 --cap 6 --out /Volumes/data/CodexTemp/T2/six73-check.json

两套源码分别位于 [128 位并集实现](../059-support/wide-engines/)及 [固定数组实现](../060-support/compact-independent/)。
