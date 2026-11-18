# R81 独立二位 bent 张量与商空间交集

[连通性证明](../../proofs/081-bent-deletions.md)说明非字符布尔谱乘独立二位满谱因子后，任意普通频率核外再删至多两点，余谱非空时的齐次拟阵仍连通。三次删除可在二维普通核断开。

| 检验 | 实际输入 | 配对数 |
| --- | --- | --- |
| Mac mini M6 CPU，Python 3.15.0 | 底谱 1–3 变量全部 276 张真值表，每张乘独立二位 bent 后枚举全部核 | 96,880，其中 90,296 个非字符配对 |
| MacBook M5 Pro CPU，Python 3.9.6 | 底谱 4、5 变量的 64、4 张固定种子表，枚举全部核 | 297,648，其中 294,823 个非字符配对 |
| Mac mini M6 CPU，Python 3.15.0 | 两个二位 bent 的张量，全部核和核外至多三删点 | 19,487，首次断连为二维核加三删点 |
| Mac mini M6 与 MacBook M5 Pro Metal | 各独立重算上述小维全部核和删点检查 | 每台 96,880 加 19,487 个配对，逐表摘要一致 |

对所有测试非字符输入均未发现未删点断连。字符底谱的首个断连维数为零或一，分别对应零频和非零频字符。MacBook 样本种子为 `810828+n`，其中 `n` 是底谱环境维数。`quotient-overlap-mini.json` 给出两变量及四变量两个真实布尔限制。在四变量例中，删除一点后外谱分为两个齐次直和块，但查询该点使两个块在真实切片上出现同一个非零商频率。该结论另用在实际限制真值表上的 Walsh 求和复核。

在归档目录内复现时沿用 `../round080/check_components.py` 的子空间枚举及基本回路判定，以及 `../round077/check_odd_walsh_metal.py` 的 Metal 运行后端。MacBook 运算经已认证 SSH 标准输入传入，输出写回 Mac mini。

```sh
PYTHONPATH=../round080 python3 -B bent_tensor.py exhaustive > bent-tensor-mini.json
PYTHONPATH=../round080 python3 -B bent_tensor.py deletions > bent-tensor-deletions-mini.json
PYTHONPATH=../round080 python3 -B quotient_overlap.py > quotient-overlap-mini.json
set -o pipefail
PYTHONPATH=../round080 python3 -B bent_tensor.py bundle ../round080/check_components.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes macbook-math /usr/bin/python3 - sample > bent-tensor-macbook.json
python3 -B bent_tensor_metal.py ../round080/check_components.py bent_tensor.py ../round077/check_odd_walsh_metal.py > bent-metal-mini.json
python3 -B bent_tensor_metal.py --bundle ../round080/check_components.py bent_tensor.py ../round077/check_odd_walsh_metal.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes macbook-math /usr/bin/python3 - > bent-metal-macbook.json
```

Mac mini 的两组 CPU 检验分别耗时 0.036 与 0.077 秒，平均占用一个核心的 98.5% 和 97.2%。MacBook CPU 检验耗时 39.853 秒，平均占用一个核心的 85.9%。Mac mini 的实际 Metal 命令区间合计 22.588 毫秒，占提交和等待时间的 90.38%。MacBook 的实际 Metal 命令区间合计 25.817 毫秒，占提交和等待时间的 89.89%。早期连接曾在服务器协议识别阶段超时。加长 SSH 握手超时后，MacBook Metal 成功运行并同步返回。

## 共同查询与选择构造

本轮 Python 3.15.0 的补充核验可从本目录执行。全部判定为二元消元、整数或有理数比较。

| 命令 | 结果 |
|---|---|
| `python3 -B check_spectral_separators.py` | spectral-separators-results.json，三、四变量分别有 168 和 29,400 个真二分隔实例 |
| `python3 -B check_overlap_cut.py` | overlap-cut-results.json，三、四变量分别穷举 1,792 与 983,040 个真值表和查询配对 |
| `python3 -B verify_interface.py` | interface.json，1,908 个实际父谱和 6,120 个限制 |
| `python3 -B verify_connectivity.py` | selector-connectivity.json，48,990 个直积与核配对以及 34 个选择树父函数 |
| `python3 -B check_aperiodic.py` | aperiodic.json，无周期选择器及各核检验 |
| `python3 -B check_near_separation.py` | near-separation.txt，6,561 张三值表与六种坐标分组，另检验门控及稀有块 |
| `python3 -B cover_differences.py` | cover-differences.json，三变量全部 30 种实际谱模长向量、1 至 4 个子空间的全部无序重复组合 |
| `python3 -B rare_gate.py mini 12` 和 MacBook 上的 `python3.15 -B rare_gate.py book 15` | rare-mini.json、rare-macbook.json，合计穷举 65,788 个三、四变量非恒种子的稀有地址构造，以及 64,516 个三位非恒种子对的乘积符号类 |

穷举采用确定性顺序。抽样用的种子见各结果文件。并行进程数及累积 CPU 时间见 rare-*.json。
