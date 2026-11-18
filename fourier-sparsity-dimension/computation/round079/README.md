# R79 仿射频层的整数核验

令 $q=|S\cap P|$，$\ell=q-\operatorname{rank}_{\mathbb F_2}\{(1,\xi):\xi\in S\cap P\}$，$t=\operatorname{codim}P$。空交集的齐次列秩取零。每个平坦均通过子空间及其陪集枚举，并按高斯二项系数校验完备性。CPU 使用整数 Walsh 蝶形变换，Metal 使用逐字符直接求和和自行实现的二元消元。

| 运行 | 输入 | 核验次数或交付物 |
| --- | --- | --- |
| Mac mini CPU | 1–4 变量全部 65,812 张真值表 | 20,132,796 个真值表和平坦的配对 |
| MacBook CPU | 5–7 变量独立固定种子样本共 548 张表 | 3,651,244 个配对 |
| Mac mini 与 MacBook Metal | 各自重算 4 变量完整枚举及 5 变量 512 张样本 | 各 21,374,464 个配对，逐表校验摘要相同 |
| Mac mini 精确回路 | 六点回路及其十个零点的全部 1,024 种填充 | 完整 Walsh 系数、支撑维数及稀疏度 |
| Mac mini 锐例 | 余维 0 至 3 的八组回路组合 | 真值填充和完整整数 Walsh 核验，最大环境维数 19 |

在归档目录内复现时，Metal 后端沿用 `../round077/check_odd_walsh_metal.py`，其 SHA-256 为 `5ef8e0baf5a142512614cbaea22511927f8262ce0cf7155e1095a363ee8a4f9a`。MacBook 的 Python 通过 SSH 标准输入运行，输出写回 Mac mini。CPU 样本种子为 `790828+n`。

```sh
python3 -B check_cut.py exhaustive > cutoff-mini.json
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes macbook-math /usr/bin/python3 - sample < check_cut.py > cutoff-macbook.json
python3 -B check_cut.py circuit > circuit-mini.json
python3 -B check_sharp_families.py > sharp-mini.json
python3 -B cut_metal.py check_cut.py ../round077/check_odd_walsh_metal.py > cut-metal-mini.json
set -o pipefail
python3 -B cut_metal.py --bundle check_cut.py ../round077/check_odd_walsh_metal.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes macbook-math /usr/bin/python3 - > cut-metal-macbook.json
```

CPU 实际运行 0.537 秒和 10.740 秒，平均占用一个核的 98.9% 和 99.9%。M6 和 M5 Pro 上的实际 Metal 命令区间累计 14.149 毫秒和 7.655 毫秒，分别占提交和同步等待区间的 79.0% 和 84.3%。GPU 时间取自 `MTLCommandBuffer`。两个主机上本次启动的运行进程均已退出。

完整随机基递推由 all_directions.py 核验。本机使用 --workers 12，MacBook 使用 --workers 15，分别生成 all-directions-mini.json 与 all-directions-macbook.json。遍历每个四变量真值表的 15 个非零方向，忽略切出空外带的方向。通过逐层最密平坦计算最优随机线性基成本，用 1 至 16 的最小公倍数作公分母；逐节点比较真实限制子函数，包括投影数据秩为零的情形。repeated.py 对具有多个共同数据标签的六回路补全按最先出现顺序的全部前驱集合和局部基作有理数求和。denominator.py 对全部四变量表按支撑和频率保留最大约分分母，逐个仿射层核对指数零度界。后两程序另以 SSH 标准输入在 MacBook Python 3.15.0 重放，结果哈希见 macbook-replay.json。

从本目录使用 python3 -B 运行 check_flat_rigidity.py、check_six.py、verify_flat.py、verify_cut.py，可分别复核任意零度取等分解、六回路补全、平坦父谱分类与纪录证明使用的小维谱形。相应结果为 flat-rigidity-results.json、six.json、flat-parent.json、verify_cut.txt。来源及无界参数证明见本轮报告。
