# 完整周期水平集的准确计算

仅需 Python 标准库（Python ≥3.9）。所有原始坐标和质量使用有理字符串；浮点输出只用于显示。`periodic_halo.py` 计算严格内界和完整外界；`independent_periodic_halo.py` 独立计算全部盒质量及矩形并面积。

从本课题目录运行，输出放 data 临时目录：

```sh
python3 -I -B computation/periodic_halo.py computation/inputs/legacy_r004_transfer_34_27.json --bits 40 --tail-radius 5 --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round001/primary.json --rectangles-out /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round001/rectangles.json
python3 -I -B computation/independent_periodic_halo.py --input computation/inputs/legacy_r004_transfer_34_27.json --bits 40 --radius-cap 5 --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round001/independent.json --compare-primary /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round001/primary.json --compare-rectangles /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round001/rectangles.json
python3 -I -B computation/tests/test_periodic_halo.py
python3 -I -B computation/verify_saved_results.py
wolframscript -file computation/symbolic_checks.wl
```

运行前确认输出目录已创建、data 已挂载，并把 TMPDIR/TMP/TEMP 及缓存指向 data。Mathematica 脚本结束显式退出。

前两条命令重新枚举几何；`verify_saved_results.py` 只核对归档结果之间的准确算术和文件身份，不重算几何。原始矩形大数组可重建，不进入公开提交。软件版本和实际参数分别记录在结果 JSON 中；不同算法使用不同 Python 版本的事实按实际保留。

一般 L¹ 下界还依赖 [完整解析传递](../proofs/periodic-complete-halo.md)，不能只凭程序输出宣称。结果中的周期比值外界不适用于全部输入，也不是 C₂ 的一般上界。

## 第2轮搜索与最终候选

`scan_perturbations.py` 和 `scan_position_grid.py` 保留实际运行版本的字节，原始默认路径对应本机实验。跨工作目录复现时使用 `replay_search.py` 显式提供路径，它不改变候选生成规则：

```sh
python3 -B computation/replay_search.py directional --input computation/inputs/legacy_r004_transfer_34_27.json --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay002/directional
python3 -B computation/replay_search.py position --input computation/inputs/r002_directional_0302.json --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay002/position
python3 -I -B computation/periodic_halo.py computation/inputs/r002_final_0475.json --bits 40 --tail-radius 5 --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay002/final-primary.json
python3 -I -B computation/independent_periodic_halo.py --input computation/inputs/r002_final_0475.json --bits 40 --radius-cap 5 --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay002/final-independent.json
python3 -I -B computation/verify_round002.py
```

完整搜索使用12个进程。只复核下界时，运行最终候选的两个算核和解析传递即可；无需重新搜索6305个点。归档的搜索账本保留各批次规模、耗时、代码/输入/任务清单哈希，候选发现的全部可再生记录留在 CodexTemp。

## 第3轮周期与联合提案

最终证书的最小核对为：

```sh
python3 -I -B computation/periodic_halo.py computation/inputs/r003_final_0546.json --bits 40 --tail-radius 5 --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay003/final-primary.json
python3 -I -B computation/independent_periodic_halo.py --input computation/inputs/r003_final_0546.json --bits 40 --radius-cap 5 --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay003/final-independent.json
python3 -I -B computation/verify_periodic_result.py --input computation/inputs/r003_final_0546.json --primary computation/results/round003/primary_final.json --independent computation/results/round003/independent_final.json --baseline-primary computation/results/round002/primary_final.json
```

最后一条只核保存算术与文件身份；前两条重算几何。搜索复现使用：

```sh
python3 -B computation/replay_search.py directional --input computation/inputs/r002_final_0475.json --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay003/directional
python3 -B computation/scan_period_grid.py --input computation/inputs/r002_final_0475.json --out /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay003/period-grid --workers 12
python3 -B computation/build_lp_input.py --results /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay003/directional/results.json --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay003/matlab-input.json
```

MATLAB 中 `addpath('computation'); propose_joint(inputPath,outputPath);` 生成数值提案，用后 `quit`。随后 `scan_joint_proposals.py --input computation/inputs/r002_final_0475.json --lp-input <提案输入> --lp-output <MATLAB输出> --out <data临时目录> --workers 12` 把联合方向有理化并完整认证。复现候选发现可直接使用已归档的 MATLAB 输入输出；不同求解器版本可能选择不同的等价 LP 解，不影响已给定准确最终输入的证明。

## 第4轮共同坐标组

`python3 -B computation/scan_coordinate_groups.py --input computation/inputs/r003_final_0546.json --out /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay004/group-search --workers 12` 重现最多八批的候选生成与准确搜索。输入原子的相等坐标用准确分数分组；当前批包含基线，后续批仅在严格区间改善时继续。

最终认证把上节两套面积命令的输入换成 `computation/inputs/r004_final_0504.json`。保存结果检查为：

```sh
python3 -I -B computation/verify_periodic_result.py --input computation/inputs/r004_final_0504.json --primary computation/results/round004/primary_final.json --independent computation/results/round004/independent_final.json --baseline-primary computation/results/round003/primary_final.json
```
