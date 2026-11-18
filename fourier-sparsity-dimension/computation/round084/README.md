# R84 非满谱的二次平移与反复提升

[解析说明](../../proofs/084-fixed-tail.md)分开处理固定零平坦的迭代、同一商的任意换基，以及每层改选零平坦的有限试验。固定零平坦时，公共直线提升的支撑与秩有精确闭式。独立地址二位块的有序迭代有权重生成函数公式。允许任意换基时，反复使用独立块仍不能产生无界比值。

CPU 在真实布尔谱中选取仿射基。对当前零集，每个贪心步骤遍历允许的全部仿射超平面来最大化仍保留的标签数，并挑出商上独立的真实谱点。用这组标签的对偶坐标构造公共直线或独立块二次相位。每一步计算是精确的整数检验。跨多个贪心步骤的结果只针对这一选择规则。自适应示例达到显式维数或支撑上限时，保留上限和后继候选记录。

| 运行 | 实际核验 |
| --- | --- |
| MacBook M5 Pro，15 个 CPU 工作进程 | 5–8 变量 16,384 张表 |
| Mac mini M6，12 个 CPU 工作进程 | 独立种子 5–8 变量 8,192 张表 |
| Mac mini 结构化多层 | 八个种子，196 个实际稀疏有符号提升，其中 99 个还独立重算完整布尔真值谱 |
| Mac mini 高秩定向试验 | 9,216 张额外表，并对实际极化秩四的明确反例重算 19 位完整真值变换 |
| Mac mini 换基试验 | 93,639 个不含零的仿射四点平面，66 个混合换基提升，以及三个第二层达到十倍活跃点数的见证 |
| M6 与 M5 Pro Metal | 各独立重算 56 次实际施加二次平移后的谱，共各 823,040 个整数系数，逐谱 SHA-256 相同 |

随机主扫描中共有 14,961 张内在非满谱表。公共直线的贪心选择在 4,584 张表上提升势，其中 472 张可选秩至少二的商。独立二位块在 4,586 张表上提升势，其中 193 张可选秩至少三的商。这些计数中的商秩与每个标签的实际极化秩分别记录。额外的 13 位构造在极化秩四下实现 $P^2=25/32\to147/160$。五位点翻号加罕见枝的稀疏例在第二层换基后实现 $P^2=27/22\to121/72\to289/144$，仍落在长期有界的模型内。

MacBook CPU 任务耗时 5.170 秒，平均使用 7.30 核。Mac mini CPU 任务耗时 1.967 秒，平均使用 5.29 核。两台 Metal 主批次实际 GPU 命令区间为 M6 43.666 毫秒、M5 Pro 31.101 毫秒。补充批次分别为 26.816 与 9.247 毫秒。主批次用对真实平移表的直接整数 Walsh 求和。补充批次的 19 位例用真种子的整数 Walsh 系数乘独立枚举的地址字符和，且已另做实际布尔全表变换。两台机器的任务进程均已退出。

Metal 运行后端沿用 `../round077/check_odd_walsh_metal.py`，SHA-256 为 `5ef8e0baf5a142512614cbaea22511927f8262ce0cf7155e1095a363ee8a4f9a`。MacBook 源码通过已认证 SSH 标准输入运行，输出写回 Mac mini。

```sh
python3 -B phase_lifts.py fixtures > fixtures-mini.json
python3 -B phase_lifts.py samples --workers 12 --samples 2048 --domain 0 > sample-mini.json
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - samples --workers 15 --samples 4096 --domain 1 < phase_lifts.py > sample-macbook.json
python3 -B high_rank.py > high-rank-mini.json
python3 -B check_rebase.py > rebase-mini.json
python3 -B phase_metal.py --witness-file high-rank-mini.json phase_lifts.py ../round077/check_odd_walsh_metal.py > phase-metal-mini.json
python3 -B phase_metal.py --extra phase_lifts.py ../round077/check_odd_walsh_metal.py > extras-metal-mini.json
set -o pipefail
python3 -B phase_metal.py --bundle --witness-file high-rank-mini.json phase_lifts.py ../round077/check_odd_walsh_metal.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - > phase-metal-macbook.json
python3 -B phase_metal.py --bundle --extra phase_lifts.py ../round077/check_odd_walsh_metal.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - > extras-metal-macbook.json
python3 -B compare_metal.py > comparison.json
```

## 外带质量、阈值及递归核验

从本目录运行以下命令，Python 3.15.0。外带质量检验穷举每个指定种子和父谱的全部齐次线性子空间，以整数质量交叉相乘比较 L。两机使用不同种子，共 2,134 次提升和 22,169,307 次子空间检查。一般二次矩阵测试另用直接 Walsh 变换重建各带。任意符号平移的仿射基测试在两机共检验 2,592 个实例。

~~~sh
python3 -B cut_mass.py 12 840240
ssh -o BatchMode=yes macbook-math /opt/homebrew/bin/python3.15 - 15 840340 < cut_mass.py
python3 -B phase_transfer.py 12 840010
ssh -o BatchMode=yes macbook-math /opt/homebrew/bin/python3.15 - 15 840100 < phase_transfer.py
python3 -B check_general_cut.py
~~~

依次对照 cut-mini.json、cut-macbook.json、phase-mini.json、phase-macbook.json 和 general-cut.json。其余单机整数验证均可用 python3 -B 加表中脚本名直接复现。

| 脚本 | 输出对照 |
|---|---|
| verify_threshold.py | threshold.json |
| verify_matching.py | matching.txt |
| verify_nonfull.py | nonfull.json |
| verify_balanced.py | balanced_verification.json |
| verify_signed_masks.py | signed-masks.json |
| check_recursive.py | recursive-results.json |
| check_unguarded.py | unguarded-results.json |

四位种子计算枚举实际 Walsh 系数并遍历频率标签的所有仿射平坦，包括空集。含全谱的最小平坦给齐次秩。对每个给定平坦维数取外质量最小值，便得到 L 的精确最大项；支撑张成之外的平坦总能与支撑张成取交，故不会高估 L。两段时间统计各自独立。

~~~sh
python3 -B four_bit_seeds.py 12 0 24576
ssh -o BatchMode=yes macbook-math /opt/homebrew/bin/python3.15 - 15 24576 65536 < four_bit_seeds.py
python3 -B verify_cut_refinement.py
~~~

输出依次见 four-bit-mini.json、four-bit-macbook.json 和 cut_refinement_verification.json。全部 65,536 张表中，恰有 512 张取最大值 8/3。
