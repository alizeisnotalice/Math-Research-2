# 五体 0.83732 外围定位和下一批门槛

将圆盘固定在原点，按三边形与五边形的二面体对称统一相对角域。原树范围内，任何面积不超过 \(0.83732\) 的五体放置必落入 [保留表](../computation/tomar-083732-retained.json) 中的 **705 个原闭叶**。表键为十进制分片编号，值为从零开始的整个原二叉树的前序节点位置。节点位置直接索引该分片的轴代码数组。候选均为原叶，分片号与节点位置共同识别每个候选。

第107、113轮、原 \(0.83731\) 四层分数及五体全域下端按 [连续门槛索引](tomar-multi-target-scores.md) 合并。面积整数单位为 \(2^{-57}\)，严格通过 \(0.83732\) 的最小整数下端为 \(120670529279675639\)。初始已排除 6,490,030 个规范原叶，其余 220,284 个进入四层。原根之外四条径向割线给出的下端严格超过 \(0.85998\)。

新四层采用能跨过 \(0.83733\) 的外截断半径。四个半径的分子为 \(1202425916136,\ 1170563102622,\ 1056804586249,\ 1003134989732\)，公共分母为 \(2^{44}\)。每次通过的整数扇和圆弓都记录整个原闭叶的下端。对受坐标截断影响的原盒，在截断内材料分数及裁掉区域的径向下端间取最小值。四层运行给出以下筛除数。

| 层 | 输入原叶 | 排除原叶 | 剩余原叶 | 节点判断 | 新整叶整数越过 0.83733 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 220,284 | 134,132 | 86,152 | 17,052,674 | 97,606 |
| 2 | 86,152 | 75,383 | 10,769 | 28,748,782 | 17,277 |
| 3 | 10,769 | 8,555 | 2,214 | 38,829,161 | 0 |
| 4 | 2,214 | 1,509 | 705 | 50,999,689 | 0 |

最终保留表有 31 个分片，与旧 \(0.83731\) 的 614 叶共有 614 个；此外有 91 个新叶。这个清单可用于任意同一第六输入。将新完整下端并入连续整数索引后，\(0.83733,\ 0.83734,\ 0.83735\) 的剩余原叶分别为 228,120、347,289、460,860。这里完成的是五体定位，第六输入的面积终端另行核验。

## 复算

MAIN 指向公开计算文件目录。TREE 为 Tomar 原 cert11d_08366。TABLE 为第113轮已认证的32等分材料表。BASE 为由旧包重新构造的 report.json 和 exact-floor 所在目录。WORK、CACHE 分别为新工作目录与 Clang 模块缓存。用相同的 \(0.83733\) 径向半径生成两个起始列表。

    python -B prepare_tomar_multi_target_layer.py BASE/report.json BASE/exact-floor \
      MAIN/five-body-radial-affine.json WORK/seeds --target .83732 --radial-target .83733 \
      > WORK/seed.json
    python -B prepare_tomar_multi_target_layer.py BASE/report.json BASE/exact-floor \
      MAIN/five-body-radial-affine.json WORK/future --target .83733 --radial-target .83733 \
      > WORK/future.json
    python -B replay_tomar_083731_scored.py MAIN TREE TABLE WORK/seeds WORK/seed.json \
      WORK/stages CACHE --workers 3
    python -B audit_tomar_083731_score_layers.py WORK/stages \
      WORK/seeds WORK/seed.json WORK/future WORK/future.json > WORK/audit.json
    python -B merge_tomar_multi_target_layer.py BASE/report.json BASE/exact-floor WORK/stages \
      WORK/seed.json WORK/audit.json WORK/floors > WORK/index.json

公开档案共 6,425,185 字节、8,194 个成员，SHA-256 为 01a2cf6ab2aa828f6aebb94667f100b9d37ad105aed8fa859bbc8f898fe20b16。可将公开档案解到 WORK/stages 以直接审计四层结果并更新整数索引。四层入口的 --parts 240 可重做一个含真实余叶的分片。在第二个已装 Python 环境独立重放了第240片，四层的逐片整数数组和配置回执字节一致。正式运行使用 macOS arm64、Python 3.13.18、NumPy 2.5.5、python-flint 0.9.0 及 Apple Clang 21.0.0。C++ 用 -ftrapv 编译。

[四层输出档案](../computation/tomar-scored-083732.zip) · [四层汇总](../computation/tomar-083732-stage-ledger.json) · [更新后的门槛索引](../computation/tomar-multi-target-after-083732.json)
