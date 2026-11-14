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

## 第5轮有限源对偶

以下复核使用标准库，无需重跑 LP 搜索：

```sh
python3 -B computation/general_dual_obstructions.py
python3 -B computation/general_dual_mass_probe.py --verify computation/results/general-dual/mass-probe.json
python3 -B computation/general_dual_saved_certificate_check.py --certificates computation/results/general-dual/mass-probe.json --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay005/independent.json
python3 -B computation/general_dual_audit_checks.py --module computation/general_dual_mass_probe.py --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay005/small.json
```

`general_dual_mass_probe.py` 不带参数时用本机 SciPy/HiGHS 重做确定种子的100例搜索，结果打印到标准输出。原探针使用 Python 3.13.18、SciPy 1.18.1；严格性由保存的准确证书和无求解器的两个检查承担。`general_dual_saved_certificate_check.py` 独立证明本批正容量饱和；它的输入契约是本批100份证书，不是任意未知容量问题的通用求解器。

## 第7轮全参数族

只复验整个参数族上界，无需重跑搜索：

```sh
python3 -B computation/independent_vertical_bound_check.py --input computation/inputs/r003_final_0546.json --certificate computation/results/round007/bb-live.json --round004 computation/results/round004/independent_final.json --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay007/independent-bb.json
```

重现分支搜索：

```sh
python3 -B computation/vertical_branch_bound.py computation/inputs/r003_final_0546.json --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay007/bb.json --target-gap 1/100000000 --seconds 180 --checkpoint-seconds 30 --cache /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay007/family.json
```

运行前创建输出父目录并确认data已挂载。墙钟上限随本机负载影响是否达到误差目标；即使时间停止，完整叶证书仍有效，实际完成状态须读 `status` 与准确缺口。加 `--resume` 从同版本代码、同输入和同目标的检查点继续。缓存不承担独立正确性证明；独立检查器从准确输入重建几何。`r007_best_vertical.json` 是该族最佳采样输入，其完整面积双算核结果位于 `results/round007/best-*.json`。

## 第8轮联合组方向

```sh
python3 -B computation/scan_joint_groups.py --input computation/inputs/r006_final_0015.json --mass-summary computation/results/round006/epoch7_summary.json --coordinate-summary computation/results/round006/epoch5_summary.json --out /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay008/joint-groups --workers 12
python3 -I -B computation/verify_periodic_result.py --input computation/inputs/r008_final_0611.json --primary computation/results/round008/primary_final.json --independent computation/results/round008/independent_final.json --baseline-primary computation/results/round006/primary_final.json
```

第一条重做3267候选搜索，第二条仅检查已保存准确算术和身份。重算最终几何时，把两套周期面积命令的输入改为 `r008_final_0611.json`，仍用 bits=40、半径上限5。方向摘要中的坐标字符串是历史标签；实际变换按成员编号作用于固定R6基准，详见本轮独立审查。

## 第9轮核族成本地板

选择已有的 Python 3.12 / python-flint 0.9.0 环境为 `PY_ARB`。下列独立最小复核不调用 LP、求逆或原 evaluator：

```sh
python3 -B computation/unpack_kernel_preconditioner.py computation/results/kernel-interface/floor-preconditioner.byteshuffle8.gz /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-interface/replay-R.bin
"$PY_ARB" -B computation/kernel_floor_neumann_review.py computation/results/kernel-interface/floor-certified.json computation/inputs/kernel-audit/candidate.json /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-interface/replay-R.bin /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/kernel-interface/replay-independent.json
```

运行前确认data挂载并创建输出目录。解码器校验压缩和原始字节两项哈希；矩阵R本身不作为可信逆矩阵，而由向外矩阵乘法验证残差。活动点必要性与完整变量列检查已包含在独立检查器中。主认证和可再生发现步骤见 [核族地板证明](../proofs/kernel-family-floor.md)；`kernel_floor_scope_check.py`另记录对全部10657行发现菜单的审计，其大菜单可再生但不属于最终证明依赖。

## 第10轮102原子超胞元

用两套原周期算核对 `computation/inputs/r010_final_0019.json` 重新计算，bits=40、半径上限8。准确保存结果检查：

```sh
python3 -I -B computation/verify_periodic_result.py --input computation/inputs/r010_final_0019.json --primary computation/results/round010/primary_final.json --independent computation/results/round010/independent_final.json --baseline-primary computation/results/round008/primary_final.json
```

搜索实现 `scan_supercells.py`、`scan_supercells_all.py`、`scan_supercells_cross.py` 保留实际字节与参数接口。搜索账本区分固定候选、缓存复用和重复夹逼；候选发现的全清单可再生，最终下界不依赖重跑它们。有限传递须用新周期和R=8所得q=4。

## 第11轮输运尾界

```sh
python3 -B computation/check_supercell_transport_tail.py --input computation/inputs/r008_final_0611.json
python3 -B computation/periodic_halo_transport.py computation/inputs/r011_transport_probe.json --base computation/inputs/r008_final_0611.json --copy-factors 1 2 --delta 1/32 --bits 40 --tail-radius 5 --output /Volumes/data/CodexTemp/T2/planar-centered-square-weak11/replay011/transport.json
```

这里1/32是合法的宽位移上限；归档实际探针使用其准确最大位移。完整几何交叉复验用原两套周期算核对同一探针、bits=40、半径上限8运行。测试脚本 `test_periodic_halo_transport.py` 的参数为含三份 `probe_*.json` 的临时结果目录，会把测试输出写入该目录；先复制相应归档JSON到data临时目录，不覆盖旧证书。独立小例及102标签复验保存于本轮 `independent-audit-checks.json`，原执行路径按实际保留。
