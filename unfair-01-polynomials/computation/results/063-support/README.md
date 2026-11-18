# 第 63 轮严格支撑范围

| 每侧最大点数 | 总次数 | 原始支撑对 | 约化支撑对 | 初筛未判 |
| ---: | ---: | ---: | ---: | ---: |
| 6 | 74 | 67,032,689,729 | 32,356,650,775 | 0 |
| 6 | 75 | 73,446,091,722 | 36,722,216,043 | 0 |
| 9 | 43 | 179,851,778,156 | 89,926,203,304 | 157 |
| 9 | 44 | 282,196,657,056 | 131,251,510,634 | 217 |

六点两套本机实现共 11,830 个分片逐条相符。MacBook 两套九点实现共 39,981 个分片逐条相符。九点的 374 个筛后余例由 368 个直接对角证书和 6 个有理恒等式排除。主进程从归档路径逐项验证，零未判。两机的已完成主管和工作进程均已退出。

    python3 -B ../056-support-compact/compare_complete.py six74-75/production six74-75/independent --totals 74 75 --cap 6 --out /Volumes/data/CodexTemp/T2/a30-six75-new.json
    python3 -B ../056-support-compact/compare_complete.py nine43-44/production nine43-44/independent --totals 43 44 --cap 9 --out /Volumes/data/CodexTemp/T2/a30-nine44-new.json
    python3 -B nine43-44/verify_residual.py

精确枚举驱动、规范代表的数学依据、对应的源路径见[第 62 轮覆盖说明](../062-support/README.md)。
