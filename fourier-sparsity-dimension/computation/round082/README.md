# R82 真实一位切片的精确搜索

每张真值表的所有非零环境方向都进入计算。方向在父谱仿射差空间内时，单独记入有效降维统计。逐方向核对父子支撑恒等式、R81 重合预算、代入真实父秩的放宽预算，以及两个子函数与父函数的 $P^2$ 比较。全部判断使用整数。随机候选按七类生成，包括均匀表、奇重量表、少数点翻转、二至四次布尔多项式，以及两支独立的选择器。

| 实际运行 | 输入 | 父表与方向配对 |
| --- | --- | --- |
| MacBook M5 Pro CPU，Python 3.9.6，15 个工作进程 | 5–8 变量，28,672 张表 | 2,658,304 |
| Mac mini M6 CPU，Python 3.15.0，12 个工作进程 | 5–8 变量，15,360 张表 | 1,360,896 |
| M6 与 M5 Pro Metal | 各自直接计算 2,048 张固定种子随机表及 18 张结构表 | 各 245,998，逐维摘要一致 |
| M6 与 M5 Pro Metal | 各自重算三个六变量反例的全部方向 | 各 189，摘要一致 |
| Mac mini 独立真切片复核 | 五个见证的真值表直接限制及变换 | 886 个实际子真值表 |

随机计算在六变量三次多项式中有 18 次出现全方向严格下降，MacBook 与 Mac mini 分别为 12 与 6 次。其余抽样维数未见相同现象。所有这 18 次也使 R81 充分判据及使用真实父秩的放宽判据逐方向失败。[六变量证明](../../proofs/082-six-bit.md)给出一个显式例，并解析计算 63 个方向的三个类别。两个独立随机任务保留的六变量见证也具有相同的精确切割统计。父 $P^2=4$，所有子函数的最大 $P^2=49/13$。这里仅判定一位限制。

CPU 并行任务根据批次设种子，主种子为 `820828`，归档的 summary.json 保留分配域、计数与运行元数据；命令可重新生成完整批次输出。Metal 独立实现直接 Walsh 求和与 GF(2) 消元，同 CPU 蝶形变换所得逐表摘要对照。单独的实际子真值表变换逐系数核对发现的例子。Metal 后端沿用 `../round077/check_odd_walsh_metal.py`，SHA-256 为 `5ef8e0baf5a142512614cbaea22511927f8262ce0cf7155e1095a363ee8a4f9a`。

```sh
python3 -B scan_cut.py structured > structured-mini.json
python3 -B scan_cut.py parallel --workers 12 --sizes 6144 4096 3072 2048 --domain 0 > parallel-mini.json
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=30 -o IPQoS=none macbook-math /usr/bin/python3 - parallel --workers 15 --sizes 10240 8192 6144 4096 --domain 1 < scan_cut.py > parallel-macbook.json
python3 -B verify_six_witnesses.py > witnesses-mini.json
python3 -B cut_metal.py scan_cut.py ../round077/check_odd_walsh_metal.py > cuts-metal-mini.json
python3 -B cut_metal.py --only-witnesses scan_cut.py ../round077/check_odd_walsh_metal.py > witnesses-metal-mini.json
set -o pipefail
python3 -B cut_metal.py --bundle scan_cut.py ../round077/check_odd_walsh_metal.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - > cuts-metal-macbook.json
python3 -B cut_metal.py --bundle --only-witnesses scan_cut.py ../round077/check_odd_walsh_metal.py | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - > witnesses-metal-macbook.json
```

MacBook CPU 任务耗时 13.955 秒，累计 187.970 CPU 秒，平均使用 13.47 核。Mac mini CPU 任务耗时 9.778 秒，平均使用 6.73 核。主 Metal 扫描的 GPU 命令区间累计为 M6 61.066 毫秒和 M5 Pro 55.828 毫秒，分别占提交和等待区间的 92.49% 和 92.76%。三例小批次的命令区间为 4.112 和 1.605 毫秒。命令时长通过 `MTLCommandBuffer` 读取。双机相关 Python 进程均已退出。

## 多位查询、匹配与相位族

在本目录运行下列脚本。除明列的抽样或匹配种子外，参数由脚本固定；秩、支撑、平方比值与证书裕量均用二元消元、整数或有理数判定。

| 命令 | 输出 | 核验对象 |
|---|---|---|
| `python3 -B check_joint_cut.py` 与 `python3 -B check_interface_two.py` | joint-cut-results.json、interface-two-results.json | 全部三、四变量的二维查询空间；固定二维接口的活动分配 |
| `python3 -B verify_multicut.py` | rank-loss.json | 196,608 个四变量父表与查询配置，917,504 个真实切片 |
| `python3 -B verify_columns.py` | columns.txt | 6,560 个地址列、二位改进例、72 个奇重量等模函数 |
| `python3 -B verify_four_lines.py` | four-lines.json | 两个不同的四线族及其共同二位查询 |
| `python3 -B verify_phases.py` 与 `python3 -B verify_rank_four.py` | phases.json、rank-four.json | 满谱平移分支的二次相位、极化矩阵及非满谱失效实例 |
| `python3 -B matching_cells.py 820001 12` 与 MacBook 上 `python3.15 -B matching_cells.py 820101 15` | matching-mini.json、matching-macbook.json | 九位地址上 54 个匹配；每个匹配枚举全部 174,762 个余维一、二单元 |
| `python3 -B mixed_matching.py 821012 12` 与 MacBook 上 `python3.15 -B mixed_matching.py 822012 15` | mixed-mini.json、mixed-macbook.json | 两个十二位地址匹配的全部一、二位混合查询证书 |

十二位匹配输出分别给出纯地址交维为 0、1、2 时的最小整数裕量。一维纯地址交枚举 8,190 个答案单元，二维交枚举 11,176,620 个。余下的混合终标签方程由每组合并的带宽下界统一涵盖。匹配由整型种子和 Python 3.15.0 的 shuffle 重建，输出含完整有序边表的 SHA-256。验证不展开具有 2,048 个终标签的完整函数真值表。

运行 `python3 -B verify_matching_recovery.py` 可读取上述两份匹配摘要、重建边表并核对各自的六问恢复证书。两例都在固定高六个地址位后得到一个无完整边的单元，其所有 64 个终标签的 4,096 个非零系数已逐项复算。见 matching-recovery.json。
同一脚本还验证两个实例均有三问提高 P 的单元，分别有 (26,460)、(28,456) 条完整边和单顶点。逐活动标签、逐地址频率复算支撑。结合单问及双问的全部混合查询证书，两例的最小不降势查询数恰为三。
