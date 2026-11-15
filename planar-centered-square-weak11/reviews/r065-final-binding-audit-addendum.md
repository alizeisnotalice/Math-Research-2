# R65 最终两候选绑定及失败诊断账本附录

2026-11-15，Asia/Shanghai。**最终来源、执行收据、完整保存覆盖、各向余项与精确质量账本核验通过。** 两份实际 V3 源证书均已存在并为 PASS，无修复；本次没有独立重算其全部源值或导数区间，没有运行 MATLAB、LP 或全 D 认证。

## 冻结理论与实现链

基础、诊断和分组增量三份正文保持已审查字节，最终再次核对：

- `r065-nonuniform-tensor-occupation-interface.md`：`fd78d51d9e1256e23df6d3ec7d13b7167cc2c601cead6f16cbb2751f82ebfcb5`。
- `r065-solver-diagnostic-scope.md`：`3b3a9818eeb82715dac35f7d98faf0d28ac8919d783ee4159064368cee92d369`。
- `r065-canonical-jump-group-interface.md`：`36b61d980fc6f3c6e6833fdaa2b1ecc53985362282df2bd8a9842748866bacbc`。

它们没有加入新的未审查数学章节；根代理后续实际结果写在本附录及轮次报告，不把历史“作者准备时未运行”误读成根代理现在没有运行。两份旧独审清单和其中全部自有 payload 也再次核对未变化：基础清单 `8ec325888eab6b61599a1296a740b51426feda1f0e915db0c1b5f11a39bc3b69`，V3 分组增量清单 `93601d2ad0758ea74b45af6e983360493842ae3e5626e5e6ef81cd44f8eafee2`。

两次实际 V3 均使用 MATLAB 源 `6121279c0386bc614bec6c949a1451016954693c3486dd9563e341f377a8d376`、精确检查器 `5124b452b656fb093a30f40b619bfbc31d7d3bb095385aa6b10f8e0fccc23e38` 和同一精确网格 `53cc33d77574f23e742811c568a49e5f5588cb7622f7de62e569e92cd5803f0d`。根目录副本、已审查作者源和 receipt/certificate 指向的哈希一致。

## 两组实际候选与收据

| 根代理目录 | h | MATLAB PID / exit | LP 保存迭代数 | 源叶／修复 | 浮点样点归一化 cost |
| --- | --- | --- | ---: | --- | ---: |
| `.local/compute/r065-v3-solve` | 3/40 | 77645 / 0 | 6158 | 1104 / 0 | 9.3144720059315134 |
| `.local/compute/r065-v3-h015` | 3/20 | 79678 / 0 | 7263 | 1104 / 0 | 5.7938065282479858 |

两次保存 solver exitflag 均为1。本次核对的是实际记录，并没有独立认证这些浮点 LP 的最优性。表中 cost 也不是一般 C₂ 上界；两份源证书都明确 `D_minimum_certified=false`。

候选／源证书 SHA256 分别为：

- h=3/40 候选 `86fcd67753fbc72714870ae3f1df3766fbd1e99d5c3ef730b3e5cb5478936576`，源证书 `af84f70fdac7fbe2d3df799a5be912537b350f536c686bcd3942695598de2cc0`。
- h=3/20 候选 `69c8a802d82e57eb46eb60cc8f570be6105699a385c00b4284e68e25cfa41f7f`，源证书 `0207380f6e59c1591c5019aa0afc71e1d0ca1ac09377ffaa341c8167f2164a1c`。

实际候选哈希、MATLAB receipt 的 output hash、源证书 input hash 完全相同。源 certificate 的 script hash 对应已审查检查器，source-receipt 命令的脚本／输入／输出路径正确，其 stdout 的 PASS 与叶数和完整证书一致。

MATLAB 两次 wall 约32.4216／29.2314秒；源检查器证书自身 elapsed 约19.6241／10.1943秒。source receipts 的进程79677／82071均 exit0，MATLAB receipts 记录的两个进程也 exit0。本次又实际 ps 核对这四个 PID 均不存在，没有启动或终止任何进程。作者归档目录 `v3-h0075-s2`／`v3-h015-s2` 各七份候选、源证书、双收据和诊断文件，合计14文件，与根代理原件逐字节相同。

## 独立覆盖、误差和非均匀质量核对

独立脚本 `computation/review_r065_final_binding.py` 从精确端点重建45个区间的中心和宽度，再枚举2025个物理胞元按 D4 分组，得到276个代表点。候选的代表点、两个宽度、orbit sizes 与 actual areas 均逐项吻合。

对每个代表矩形按两个实际宽度各二分，独立生成全部 `276×2²=1104` 盒。每份 accepted leaves 带胞元编号的有理端点多重集与其准确相等，所有叶 depth0，无重复漏盒。两份合计2208盒：checked=accepted，unresolved、unchecked、splits、负点数组均为空，保存 lower_fixed 全非负、repair_upper=0。

逐叶从保存的 signed curvature intervals 上端重新计算

    ceil(max(0,Uxx)*ell_x²/8)+ceil(max(0,Uyy)*ell_y²/8)，

与保存 error_fixed 完全相同。这核对了真实横纵边长的各向余项账本，**没有重新求 g 顶点值或重新验证全部导数区间**。其源数值正性来自已审查、已绑定检查器的根代理实际运行，不能写成本次独立全源重算。

将保存十进制作为精确有理数，按每个物理矩形面积独立重算

    M=lambda+Σ_physical a_orbit*width_i*width_j，
    halfmass/pi=h*M。

结果与源证书准确相同：

| h | 精确 M | 精确 halfmass/pi |
| --- | --- | --- |
| 3/40 | `506007641441308124843/12800000000000000000` | `1518022924323924374529/512000000000000000000` |
| 3/20 | `1573739431746150174497/128000000000000000000` | `4721218295238450523491/2560000000000000000000` |

矩阵诊断也对应实际新表示：V3 3543变量、529跳变等式、4416源行、6880总不等式行；未把V2的10765变量／4140跳变等式文件混入新候选。

## 三份历史诊断的边界保持

1. 原V1：MATLAB PID66280，exit1，source hash5f00c7f4...，当前进程不存在。原目录没有 solver-status；真实 exitflag／返回原因仍为 UNKNOWN，不能用后面的V2结果回填。
2. λ-only V2：PID68608，exit0，code3d696e00...，当前不存在。4416行最小值 `3.4273796974451733e−5` 为保存的浮点行诊断；网格／h／subdiv和输出哈希吻合。本次没有将其升级为精确有限行或全 D 可行证书。
3. V2 solve：PID70469，MATLAB exit1、当前不存在；保存 solver exitflag=0、iterations=20000、无返回向量／目标，属于实际迭代上限，不是数学不可行证书。矩阵10765变量／4140等式、4416源行与V2参数一致。

以上三份 receipt 的 source hashes 与对应原代码、精确网格副本均匹配。连同两份新MATLAB与两份Python进程，共七个相应 PID 的当前缺席均已实际核对。只陈述这些指定进程，不对其他用户进程作判断。

## 最终范围

本次实际 Python3.13.18，0.089356秒、exit0，仅执行记录、分片、质量和保存误差算术核验。输出 `computation/results/round065-final-binding-independent.json` 保存所有实际来源、哈希、进程状态和数量。

最终统一结论：非均匀张量数学接口与规范分组已独审；两次V3候选源PASS分别覆盖1104叶，0修复，记录与来源绑定通过；本次没有重跑所有源值。全 D 认证、改善一般 C₂ 界及精确 C₂ 仍未解决。未提交、推送或运行 MATLAB，根代理负责统一归档与发布。
