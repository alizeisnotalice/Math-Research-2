# R80 外谱拟阵连通块的整数核验

对每张布尔真值表计算整数 Walsh 变换。对每个线性子空间 $H$，在外谱 $S\setminus H$ 的齐次列 $(1,\xi)$ 上进行二元消元，以基础回路形成连通块。核对块数至多三。一位切割核对块数至多二。若恰有三块，逐块核对该块在模 $H$ 的投影只含一个标签，三个标签均非零且和为零。子空间目录按逐维高斯二项系数校验完备性。

| 实际运行 | 核验范围 | 核验数 |
| --- | --- | --- |
| Mac mini M6 CPU，Python 3.15.0 | 1–4 变量全部 65,812 张真值表和全部线性子空间 | 4,395,096 个配对，其中 79,192 个三块情形 |
| MacBook M5 Pro CPU，Python 3.9.6 | 5–7 变量固定种子抽样 292 张真值表和全部线性子空间 | 302,992 个配对 |
| M6 和 M5 Pro Metal | 各自直接重算 4 变量全部真值表、5 变量 256 张样本及所有一位切割 | 各 990,976 个配对，逐表校验摘要与 CPU 一致 |
| Mac mini 整数构造 | 2、4、6、8 变量各一例 | 三块大小依次为 $(1,1,1)$、$(4,1,1)$、$(4,4,1)$、$(4,4,4)$，全谱稀疏度最后为 76 |

CPU Walsh 变换使用整数蝶形算法，Metal 通过直接逐项求和生成支撑，并独立实现基础回路并查集。GPU 后端沿用 `../round077/check_odd_walsh_metal.py`，SHA-256 为 `5ef8e0baf5a142512614cbaea22511927f8262ce0cf7155e1095a363ee8a4f9a`。CPU 样本种子为 `800828+n`。MacBook 以 SSH 标准输入运行，结果写回 Mac mini。

```sh
python3 -B check_components.py exhaustive > components-mini.json
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes macbook-math /usr/bin/python3 - sample < check_components.py > components-macbook.json
python3 -B check_triple_fixture.py > triple-fixtures-mini.json
python3 -B components_metal.py check_components.py ../round077/check_odd_walsh_metal.py > components-metal-mini.json
set -o pipefail
python3 -B components_metal.py --bundle check_components.py ../round077/check_odd_walsh_metal.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes macbook-math /usr/bin/python3 - > components-metal-macbook.json
```

两台 CPU 分别运行 0.369 秒和 14.992 秒，平均使用单核的 99.88% 和 99.98%。`MTLCommandBuffer` 在 M6 和 M5 Pro 上记录到 6.739 和 9.261 毫秒的实际 GPU 命令区间，分别占提交和同步等待区间的 83.10% 和 76.82%。两台机器上本次启动的计算进程已退出。

## 两机穷举四维的所有符号配置

另一个独立 CPU 作业检验任意非零简单二元列上的分层切割。枚举四维空间 15 个非零点的所有非空子集。对每个支撑先固定一个贪心基，将其标记归零；模掉线性泛函只置换候选切割，所有标记均被覆盖。对每种代表、每个并列最优切割、每个密度层和每个含前一层的子空间检验半数命中与秩差容量。MacBook 以 15 个工作进程覆盖 9,565,938 种标记，Mac mini 以 12 个工作进程覆盖 4,782,968 种，合为 3¹⁵−1。总计 898,856 个符号轨道，检查 1,488,645 个并列切割。两个分片的 CPU 时间及墙钟时间分别见 exhaustive-macbook.json 与 exhaustive-mini.json。两台机器另以不同种子为每个支撑各抽一种标记，单独复验。

以下命令分别在已有这两个 Python 文件的对应主机运行，Python 3.15.0；各脚本只写标准输出。

```sh
python3 -B exhaustive_signs.py 15 book
python3 -B exhaustive_signs.py 12 mini
python3 -B signed_flags.py 15 8002
python3 -B signed_flags.py 12 8001
python3 -B check_outer_components.py
python3 -B verify_packing.py
python3 -B check_separated.py
```

后面三个脚本分别复核齐次空间中的核、成对商谱密度及完整单带剥离。对应输出依次为 outer-components-results.json、packing.json、separated.json。计算均为整数或有理数，无容差。
