# R86 Petersen 中心割与同时保秩扩展

[解析记录](../../proofs/086-petersen.md)使用十五条 Petersen 边槽和独立十四位 bent 稠密槽。相位矩阵取自顶点矩阵的二元差，稠密槽保持相位零。有效新地址维数至多 31 时，未包含整个稠密槽张成的割已被稠密质量支付。其余割化为 8,581 个图连通闭包分划。对中心相位，下界直接给出真实中心割。

| 运行 | 结果 |
| --- | --- |
| MacBook M5 Pro CPU，15 个进程 | 64,824 个横向余维至少八的真实谱矩阵族，39,163 个混合秩二和秩四，10,461 个逐遍历全部分划 |
| Mac mini M6 CPU，12 个进程 | 10,358 个合格矩阵族，其中 6,251 个混合秩二、四，1,674 个完整分划优化 |
| 两台 Metal 完整闭包割 | 各 1,298 个矩阵族乘 8,581 个分划，共各 11,138,138 个配对。逐矩阵族最小预算相同，摘要一致 |
| 两台 Metal 真实有符号谱 | 各重算 262,144 个稠密槽有符号系数及五个实际父函数的 46,080 个图边相位整数系数，摘要一致 |
| Mac mini 扩展 | 对 5,205 个额外矩阵族计算精确兼容维数。对靠近门槛及不同非局部商维数的 24 个代表做同时扩展检验 |

全部 75,182 个主筛通过可证明的分划余量下界。另有 12,135 个完全优化，共 104,130,435 个矩阵族和分划配对。混合秩中心谱的最小实际余量为二。主筛对其中 13,972 个矩阵族解出兼容线性空间，1,165 个允许至少一轴非局部保秩扩展。最接近边界的已扩展示例在十维目标恰取等。一个混合秩例的非局部商为二维，但六条回路上的配对综合征非零，最多只能同时增加一轴。较大商的扩展示例和已穷举最大各向同性的结果分别标记。未发现真实负余量。

MacBook 主 CPU 耗时 5.537 秒，累计 51.929 CPU 秒，工作阶段平均使用 9.649 个核心。Mac mini 主 CPU 耗时 0.894 秒，平均使用 5.978 核。完整割 Metal 命令区间为 M5 Pro 2.513 毫秒和 M6 6.508 毫秒。有符号谱另用 1.744 毫秒和 4.688 毫秒。命令时长由 `MTLCommandBuffer` 读取。计算进程在两台机器上均已退出。

Metal 后端沿用 `../round077/check_odd_walsh_metal.py`。最接近边界的见证与小于零的筛选均使用精确整数或有理数。样本固定种子。用 SSH 标准输入运行 MacBook 版本并写回输出。随机扫描中若多份矩阵等号并列，所留见证会随进程完成顺序变化。固定已归档矩阵输入可复验相同签名。

```sh
python3 -B petersen.py --workers 12 --tasks 48 --draws 256 --domain 0 > petersen-mini.json
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - --workers 15 --tasks 75 --draws 1024 --domain 1 < petersen.py > petersen-macbook.json
python3 -B simultaneous.py > simultaneous-mini.json
python3 -B signed_witnesses.py petersen.py ../round077/check_odd_walsh_metal.py petersen-macbook.json simultaneous-mini.json > truth-witnesses-mini.json
python3 -B graph_metal.py --pencil-file petersen-macbook.json petersen.py ../round077/check_odd_walsh_metal.py > graph-metal-mini.json
python3 -B signed_witnesses.py --metal petersen.py ../round077/check_odd_walsh_metal.py petersen-macbook.json simultaneous-mini.json > truth-metal-mini.json
set -o pipefail
python3 -B graph_metal.py --bundle --pencil-file petersen-macbook.json petersen.py ../round077/check_odd_walsh_metal.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - > graph-metal-macbook.json
python3 -B signed_witnesses.py --bundle petersen.py ../round077/check_odd_walsh_metal.py petersen-macbook.json simultaneous-mini.json | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - > truth-metal-macbook.json
```

## 其他精确核验

在本目录执行以下命令，各脚本将结果输出到标准输出。除 verify_cyclic_gap.py 需要 SciPy 外，其余证明核验使用 Python 标准库。

| 程序 | 归档输出 |
| --- | --- |
| verify_unrestricted.py | quadratic.json |
| check_unrestricted.py | audit-results.json |
| check_support_capacity.py | support-results.json |
| check.py | record.json |
| check-many.py | result-many.json |
| verify_cyclic_gap.py | cyclic.json |
| verify.py | verify.txt |
| support_orbits.py | support-orbits.json |

四变量割核验按全局变号配对，再把 32,768 个真值表代表分给两台机器。每项的 Walsh 变换和全部有效子空间上的余维均用整数计算。输出的 ratio 为最佳常值权重的二范数除以 P。独立的 support_orbits.py 用 Hamming 距离重算支撑，按可逆仿射变换生成轨道，逐轨用有理割乘子同时验证原问题可行性和对偶取等。

~~~sh
python3 -B support_census.py 0 8192 12 > census-mini.json
ssh -o BatchMode=yes macbook-math /opt/homebrew/bin/python3.15 -B - 8192 32768 15 < support_census.py > census-macbook.json
python3 -B support_orbits.py > support-orbits.json
~~~
