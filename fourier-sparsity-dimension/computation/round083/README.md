# R83 匹配选择器的共同一问与两问

[解析证明](../../proofs/083-finite-matching.md)讨论独立私有二位 bent 种子，既刻画 R49 高度一的七个纯地址恢复平面，也证明满足大混合性的匹配可在父比值趋于 3/2 时排除全部一问、两问。

数值核验使用内在频率差空间。按询问空间和纯地址、纯数据空间的交维，以及询问空间所含实际父谱差方向数，记录全部答案的支撑大小与齐次秩。含足够多谱差方向的空间用真实有符号谱合并逐项计算。一条活跃差方向的二维空间按第一问子谱在第二问方向上的秩变化精确归类计数。没有活跃差方向的情形，由坐标交维的高斯二项式计数给出。

| 结构 | 内在维数 | 父比值平方 | 一问最优平方 | 两问最优平方及恢复数 |
| --- | --- | --- | --- | --- |
| R49 字符叶 | 6 | 4 | 49/13 | 4，31 |
| R49 张量同一独立 bent | 8 | 25/16 | 4 | 4，10,731 |
| R49 私有无守卫 bent | 11 | 169/49 | 36/11 | 4，7 |
| R49 私有有守卫 bent | 14 | 4 | 49/13 | 4，恰好七个纯地址面 |
| 重复边方向，私有有守卫 bent | 14 | 4 | 4 | 4，8,195 |
| R49 私有有守卫三位奇重量种子 | 18 | 25/8 | 361/112 | 27/8，在经最优一问的 131,071 个平面中最优，有 84 个严格改善面 |

前五族的一、二维内在问空间均精确计数。三位奇重量种子枚举所有一维问空间和通过固定方向 66 的所有二维面，另直接核对全部七个纯地址面。前五族与奇重量的一问、固定方向两问总账为 90,599,920 个问空间出现次数。独立重构和变换 114 张实际子真值表，逐系数核对最优见证及 R49 的纯地址平面。

MacBook CPU 以 15 个工作进程完成两种有守卫 bent 的主扫描，耗时 7.778 秒，平均使用 11.91 核。Mac mini 的 12 进程主扫描耗时 1.926 秒，平均使用 4.81 核。两台 Metal 独立重算各 590,448 个有符号查询，完整类型直方图与 CPU 相同，两台逐问输出摘要一致。M6 及 M5 Pro 主 Metal 命令区间分别为 43.836 与 17.260 毫秒，奇重量单问批次另用 0.723 与 0.302 毫秒。命令时长来自 `MTLCommandBuffer`，与设备整体利用率不同。

Metal 后端使用 `../round077/check_odd_walsh_metal.py`，SHA-256 为 `5ef8e0baf5a142512614cbaea22511927f8262ce0cf7155e1095a363ee8a4f9a`。MacBook 通过已验证的 SSH 连接从标准输入运行，实际结果写回本机。双机相关进程均已退出。

```sh
python3 -B matching_queries.py --workers 12 r49:character r49:global-bent r49:private-bent > matching-mini.json
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - --workers 15 r49:private-guarded-bent repeated:private-guarded-bent < matching_queries.py > matching-macbook.json
python3 -B matching_queries.py --workers 12 --depth 1 r49:private-guarded-odd > guarded-odd-mini.json
python3 -B check_witnesses.py > actual-witnesses-mini.json
python3 -B matching_metal.py matching_queries.py ../round077/check_odd_walsh_metal.py r49:character r49:global-bent r49:private-bent r49:private-guarded-bent repeated:private-guarded-bent > matching-metal-mini.json
python3 -B matching_metal.py --depth 1 matching_queries.py ../round077/check_odd_walsh_metal.py r49:private-guarded-odd > guarded-odd-metal-mini.json
set -o pipefail
python3 -B matching_metal.py --bundle matching_queries.py ../round077/check_odd_walsh_metal.py r49:character r49:global-bent r49:private-bent r49:private-guarded-bent repeated:private-guarded-bent | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - > matching-metal-macbook.json
python3 -B matching_metal.py --bundle --depth 1 matching_queries.py ../round077/check_odd_walsh_metal.py r49:private-guarded-odd | ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=45 -o IPQoS=none macbook-math /usr/bin/python3 - > guarded-odd-metal-macbook.json
python3 -B compare_runs.py > comparison.json
```

## 渐近定理的有限检验

从本目录用 Python 3.15 运行下列脚本。除运行时长及主机信息外，结果可直接与相应归档比较。

| 脚本 | 对照记录 |
|---|---|
| verify_localization.py | localization.json |
| check_thinned.py | thinned.json |
| verify_bent_depth.py | bent-depth.txt |
| verify_all_ranks.py | quadratic-phases.json |
| check_one_three.py | one-three-results.json |
| check_sharp_loss.py | sharp-loss-results.json |
| check_zero_band.py | zero-band-results.json |

随机匹配检验在两台机器上使用互不相交的固定种子。每个样本遍历全部 131,070 个地址仿射超平面，并以整数和有理数检查嵌套单元的内边数及末端恢复。复现命令如下。

~~~sh
python3 -B matching_profile.py 12 831001
ssh -o BatchMode=yes macbook-math /opt/homebrew/bin/python3.15 - 15 831201 < matching_profile.py
~~~

相应记录为 profile-mini.json 和 profile-macbook.json。
