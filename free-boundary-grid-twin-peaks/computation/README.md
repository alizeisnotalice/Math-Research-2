# 计算复现

各轮使用 Python 标准库、C++ 整数程序或 Swift/Metal；具体依赖、整数位宽和适用尺寸以各程序及轮次证书为准，不得擅自放大。随机搜索不承担渐近证明或混合认证。

在 `/Volumes/data/CodexProjects/T2` 工作。依赖data写入前确认其已挂载。临时输出根为 `/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks`，C++共享缓存为 `/Volumes/data/CodexTemp/cache/clang`。

## 核验已归档证书

```sh
python3 -B free-boundary-grid-twin-peaks/computation/verify_results.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_insertion_n3.py
python3 -B free-boundary-grid-twin-peaks/computation/compute_global_geodesic_p.py
python3 -B free-boundary-grid-twin-peaks/computation/thin_cycle.py
```

## 重建精确计数

```sh
python3 -B free-boundary-grid-twin-peaks/computation/exact.py --max-n 4 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/exact-rerun.json
```

每个C++源用下列方式编译（把moment_dp换为cut_dp、pair_dp、geodesic_dp、merge_tail或conditional_walk）：

```sh
TMPDIR=/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang++ -std=c++17 -O3 -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang free-boundary-grid-twin-peaks/computation/moment_dp.cpp -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/moment_dp
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/moment_dp 5
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/geodesic_dp 5
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/merge_tail 5 8
```

merge_tail的第二参数t给严格尾事件M>t；将t遍历2至25即可复现完整尾部。峰对与切割扫描：

```sh
python3 -B free-boundary-grid-twin-peaks/computation/pair_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/pair_dp --workers 3 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/pairs-rerun.json
python3 -B free-boundary-grid-twin-peaks/computation/cut_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/cut_dp --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/cuts-rerun.json
```

## 随机实验

```sh
python3 -B free-boundary-grid-twin-peaks/computation/search.py --steps 200000 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/search-rerun.json
python3 -B free-boundary-grid-twin-peaks/computation/walk_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/conditional_walk --steps 5000000 --workers 3 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/walks-rerun.json
python3 -B free-boundary-grid-twin-peaks/computation/walk_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/conditional_walk --steps 5000000 --workers 3 --reflect-every 10 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/walks-reflect-rerun.json
python3 -B free-boundary-grid-twin-peaks/computation/walk_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/conditional_walk --steps 20000000 --workers 3 --reflect-every 10 --top-block-every 4 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/walks-block-rerun.json
```

种子由源码固定并逐链记录。walk前20%作烧入，之后每n²步记录一次；保留前后半段直方图和独立核验的最佳候选。提议包括固定分布的标签/位置交换、可选条带对合，以及可选含最高峰2×2块的均匀置换。块由最高峰位置选择时，接受比为该点原有可选块数/新位置可选块数；其他提议对称。所有更新均要求最终恰有两峰。均匀条件律平稳不等于不可约或混合已证。

浮点仅用于退火接受率与探索摘要；精确证书不依赖这些值。不同C++标准库实现的shuffle可能改变随机轨迹，实际编译器版本已记录。旧阶段输出保留当时二进制哈希；最终源码通过相应参数保留三阶段算法，性能计时不是复现相等项。

## 第2轮新增核验

`weighted_promotion.py` 完整枚举3×3，验证所有提升关系元、逐输出归一化、精确逆像及有理对数证书。`promotion_dp.cpp` 接收边长n（2≤n≤5），以第一触邻递推给出总提升质量的精确分子/分母；n5内存1 GiB。编译缓存沿用上节，命令和实际结果在rounds/002.md与manifest-r002.json。

## 第3轮新增核验

```sh
python3 -B free-boundary-grid-twin-peaks/computation/verify_r003.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_shortest_lift.py
python3 -B free-boundary-grid-twin-peaks/computation/beta_fiber_n3.py
python3 -B free-boundary-grid-twin-peaks/computation/geodesic_fiber.py --n 4 5 6 7 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/geodesic-fiber-rerun.json
```

按上节同样编译singleton_dp.cpp与beta_walk.cpp；前者n5内存1 GiB，命令为`singleton_dp 5`。后者新增任一峰2×2块和β观测，复现命令：

```sh
python3 -B free-boundary-grid-twin-peaks/computation/walk_scan.py --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/beta_walk --steps 20000000 --workers 3 --reflect-every 10 --top-block-every 4 --output /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/beta-walks-rerun.json
```

此处沿用驱动参数名top-block-every，但beta_walk中从两峰等概率选择。给定块B，选择概率正比于Σ_{峰v∈B}1/b(v)，b(v)是含v的2×2块数；接受比为输出/输入的此项之比。仅理想随机核的平稳性有解析依据，实际浮点/伪随机实现及有限轨迹不提供认证概率。全部最佳标号以独立峰判定器再次核验。

## 第4轮：独立集、缩图与全n构造证书

按上节C++编译方式构建`independent_peak_dp.cpp`（仅接受2≤n≤6），然后运行：

```sh
python3 -B free-boundary-grid-twin-peaks/computation/verify_independent_peak_dp.py /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/independent_peak_dp
python3 -B free-boundary-grid-twin-peaks/computation/check_selected_joint_peaks.py /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/independent_peak_dp
python3 -B free-boundary-grid-twin-peaks/computation/verify_r004.py
python3 -B free-boundary-grid-twin-peaks/computation/rectangle_root_counts.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_orbit_core.py
python3 -B free-boundary-grid-twin-peaks/computation/beta_zero_recurrence_n3.py 3
```

n6独立集表5598861项，主Agent实测最大常驻约532MiB；C++表用192位、交错和用256位，解析界见证明。selected checker以Python任意精度重算两个18点独立集的全子集。n6完整JSON来自实际C++输出，不手动删去根分布。

β有向边轨道程序需要本机GMP 6.3.0，实际编译方式（安装路径若变化应从本机查实）：

```sh
TMPDIR=/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang++ -std=c++17 -O3 -I/opt/homebrew/Cellar/gmp/6.3.0/include -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang free-boundary-grid-twin-peaks/computation/beta_zero_orbits.cpp -L/opt/homebrew/Cellar/gmp/6.3.0/lib -lgmpxx -lgmp -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/beta_zero_orbits
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/beta_zero_orbits 6
```

程序逐轨释放三张整数表，最大单轨773426项；每步检查整除、非负、反向边对称及有向边覆盖总数。拒绝n=1或n≥7。未测量该程序峰值内存，不把状态数当实际内存测量。

`rectangle_root_mod.cpp`接受行数和列数，限面积≤24；以模2^64和模(2^64−1)两遍连通顺序递推及CRT恢复，结果≤24!小于两模乘积，恢复唯一。它与Python独立集根带算法交叉。`core_broken_grid.cpp`无参数，完整遍历3686400个有限核心状态；独立Python验证器利用删边单调性复核，亦生成n5/6/12/32的完整近峰成员。

计算结束先保存必要输出，再确认所启动进程退出；不保留大表常驻REPL。n7采用新算法和独立资源预算，不能直接放宽旧程序输入限制。

## 第5轮：n7 D4整数计数与分轨认证

`n7_d4_peak_dp.cpp` 无参数，只处理n7；`d4_peak_dp_generic.cpp` 接收n=3..6。
使用前述data共享缓存，以 `clang++ -std=c++17 -O3 -pthread` 编译。n7主计数最高实测footprint约6.11GB。

```sh
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/n7_d4_peak_dp > /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/n7-rerun.json
python3 -B free-boundary-grid-twin-peaks/computation/check_selected_n7.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_d4_orbits.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_r005.py
```

选定集合脚本读取归档n7结果，使用Python任意精度重算。Burnside脚本使用独立的周期商图与逐行算法，不重复C++D4枚举。
`rectangle_independent_root.cpp` 接收行数、列数，限两者≤7且面积≤35；`rectangle_independent_root 7 5` 复现最大矩形。

β完整认证在 `certificates/r005-beta/`，保留冻结源码、顺序驱动、输入、每轨JSON/资源和必要历史对照。该源码同样以C++17、O3、pthread编译，不需GMP；每轨内部12线程。不要同时启动全部12轨，驱动保证上一轨退出后才启动下一轨。实际最大footprint约8.78GB。

```sh
python3 -B free-boundary-grid-twin-peaks/computation/certificates/r005-beta/verify_beta_zero_layered.py --orbit-dir free-boundary-grid-twin-peaks/computation/certificates/r005-beta/orbits --certificate /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/beta-check.json
python3 -B free-boundary-grid-twin-peaks/computation/certificates/r005-beta/run_n7_beta_orbits_certified.py --source free-boundary-grid-twin-peaks/computation/certificates/r005-beta/beta_zero_layered_orbit_certified.cpp --binary /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/beta_zero_layered_orbit_certified --output-dir /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/beta-rerun
```

第二条命令会实际重跑12轨，输出目录必须不存在。历史单进程的binary哈希缺口与新分轨认证严格分开；见证书目录SOURCE_VERSION.md。各进程完整退出后释放表和线程，不保留数值REPL。

## 第6轮：前二扇区与Metal搜索

`n7_top2.cpp` 无参数，沿用n7 D4表并新增前二扇区两个精确累加器；以C++17、O3、pthread编译，实际footprint约5.06GB。

```sh
python3 -B free-boundary-grid-twin-peaks/computation/audit_top2_identity.py
python3 -B free-boundary-grid-twin-peaks/computation/audit_top2_first_merge.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_support_T_n3.py
python3 -B free-boundary-grid-twin-peaks/computation/search_two_step_n4.py
python3 -B free-boundary-grid-twin-peaks/computation/n8_independent_set_census.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_r006.py
```

`top2_h_dp.cpp`接收n=3..5，以反向增长DP输出逐种子对h；`audit_top2_h.py`接受该binary路径作Python任意精度独立包含—排除核验（参数见源码）。n8脚本仅做逐行状态数资源预检，未枚举n8的C表。

Metal源码、原始种子/参数记录与完整样本在 `certificates/r006-gpu/`。用Swift6.4.2及Metal运行时编译，不要求离线metal编译器：

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-modules swiftc -O -module-cache-path /Volumes/data/CodexTemp/cache/swift-modules free-boundary-grid-twin-peaks/computation/certificates/r006-gpu/twin_peaks_gpu.swift -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/twin_peaks_gpu
```

从临时输出目录运行binary，参数为 `8 12 16 24 32` 或 `--cut-only 24 32`；该历史源码会写当前目录gpu-search-results.jsonl，重跑前必须选新目录以免覆盖旧记录。最终源码含接受次数计数，旧V1/V2记录不能冒称同一binary轨迹。历史 `gpu_seconds` 字段是提交等待墙钟，不是硬件执行秒数。
`annotate_samples.py 输入.jsonl 输出.jsonl`以CPU独立核验排列、峰数、距离、β、首次合并和界面。全209保存样本实际复核通过；AGX满载采样是系统字段，不提供混合或概率认证。完成后退出独立进程释放Metal资源。

## 第7轮：高β族归一化与完整n5两步核

```sh
python3 -B free-boundary-grid-twin-peaks/computation/verify_geodesic_seed_family.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_beta_one_n3.py
python3 -B free-boundary-grid-twin-peaks/computation/certificates/r007-gpu/high_beta_seed_audit.py
```

`two_step_completion_scan.cpp`按C++17、O3、pthread方式编译，参数为n（2..5）；n5约512MiB表、12线程。
`certificates/r007-gpu/`保存真实v1/v2源码、最终16标签、独立审计、硬件计时与版本更正。`audit_r007.py 输入.jsonl 输出.json`复核完整样本。历史源需要当前目录存在空的r007-metal.jsonl且会覆盖它，复现必须使用全新临时目录。真实种子、旧JSON错误及缺失早期样本均在README说明，不补造记录。所有程序完成后退出释放内存。

## 第8轮：独立整数反例

```sh
c++ -O3 -std=c++20 -pthread -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang free-boundary-grid-twin-peaks/computation/completion_ie_n8_largeS.cpp -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/kernel8
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/kernel8 0x1010101010107fc
A022_CERT_DIR=/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/n8-independent-python python3 -B free-boundary-grid-twin-peaks/computation/verify_largeS.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_r008.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_peak_polynomial_counterexample.py
python3 -B free-boundary-grid-twin-peaks/computation/check_extended_l_holes.py
python3 -B free-boundary-grid-twin-peaks/computation/certificates/r008-gpu/summarize.py
```

输入与两份独立整数结果在results/n8-kernel-*-r008.json。C++约375MB、Python约365MB；结果保存后进程正常退出。n7选定形态源码completion_ie_selected.cpp，参数n、S掩码（其driver保存原运行路径，迁移复现需改为本机新临时目录）。GPU入口与计时/边界更正见certificates/r008-gpu/README.md。

## 第9轮：原条件律计重与n9修正复算

```sh
python3 -B free-boundary-grid-twin-peaks/computation/original_law_tree_audit.py 3
python3 -B free-boundary-grid-twin-peaks/computation/verify_minimum_deletion.py --n 3
python3 -B free-boundary-grid-twin-peaks/computation/check_inert_reduction.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_inert_fibres.py
python3 -B free-boundary-grid-twin-peaks/computation/verify_r009.py
```

n9源码与完整generate/calculate/audit命令在certificates/r009-n9/RUN_SUMMARY.md；必须不用-NDEBUG，输出放全新data临时目录。计算约5GB RSS、外存约12GB；本轮实际修正重跑已完成，不需仅为读取结果重复运行。层检查脚本会在自身目录写同名JSON，若要保持归档只读可复制到新临时目录运行。GPU旧partial证据仅1195条，详见r009-gpu-partial/README.md。

## 第10轮

自然森林、界面和安全根：`audit_natural_forest.py`、`audit_interface_weight.py`、`audit_safe_root_forcing.py`（输出路径参数应放 CodexTemp）。全 safe 切割及 GPU/CPU 完成批次见 `certificates/r010-*`，各包附源码及复现说明。执行 `python3 -B computation/verify_r010.py`（从课题根目录）重算精确账本和数值汇总；其通过不表示原问题闭合。当前文件哈希见 `manifest-r010.json`。

## 第11轮

`audit_interface_gluing.py` 与 `audit_interface_theta_family.py` 为路径/全交错有限精确核验，输出参数请指向CodexTemp。`verify_r011_cpu.py`、`verify_r011_gpu.py` 从已归档原始批次重算矩与数值诊断。两类生产包见 `certificates/r011-cpu-variance`、`certificates/r011-gpu-sis`；实际版本、seed、资源、复现命令与限制各自记录。数学反例和归一化排除见 `proofs/interface-gluing-sharpness.md`。

## 第12轮

`audit_tau_smc_aggregates.py` 独立合并原条件律粒子重复并逐值核对GPU选父迁移；`audit_integer_stratified.py` 核验整数分层引理；`verify_layer_split_tree.py` 核验全n层分裂树的有限实例。新证书在课题根目录 `certificates/r012-first-merge/`、`certificates/r012-tau-smc/`、`certificates/r012-gpu-parent/`。均值/方差定点端点以整数分子分母为准。

第13轮：`audit_late_merge_families.py` 核验晚合并两族；`audit_r013_gpu.py` 独立审计新核来源、4800个实际seed、原条件矩与整机遥测。证书在课题根 `certificates/r013-analytic/` 和 `certificates/r013-gpu/`。

第14轮：`audit_early_interface.py` 检查可逆交错及完整方格单向阶乘例；`audit_r014_gpu.py` 以Decimal60重算n8/n9矩、一峰归一化和明确同P合并。课题根 `certificates/r014-boundary/` 保存出生及有理归一化核验，`r014-dominating/` 为n3–n5支配条件整数DP，`r014-completion/` 为n3/n4一般双向完成数和严格dyadic求和区间；各包附复现命令。小粒子大掩码核审计在 `r014-gpu/kernel-audit/`。编译二进制和重跑输出均应写CodexTemp；冻结源中的绝对临时输出路径在包README注明。


## 第15轮复现

在课题根目录运行，二进制和缓存均放data临时目录。数学参数与结果见rounds/015.md。

```sh
cd certificates/r015-completion
clang++ -O3 -std=c++20 -pthread -fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang completion_gain_compact.cpp -I/opt/homebrew/opt/gmp/include -L/opt/homebrew/opt/gmp/lib -lgmpxx -lgmp -o /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r015-completion
/Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r015-completion 5 /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r015-completion-n5.json 4
python3 verify_r015.py
python3 audit_r015_n3.py
cd ../r015-hook
python3 audit_recursive_hook_tree.py /Volumes/data/CodexTemp/T2/free-boundary-grid-twin-peaks/r015-hook.json
cd ../r015-t-interface
python3 audit_square_dstar.py
python3 verify_t_interface.py
python3 audit_bounded_parent.py
python3 t_family_exact.py
```

`root-reproduction.json`记录重新编译后n5输出除运行时间字段以外逐项完全一致。
T族分母沿用本课题已核验n3/n5/n7整数A₂；未新做n7全图计数。


## 第16轮复现

GPU完成批的版本、seed、命令及资源记录见 `certificates/r016-gpu/README.md`；`audit_r016_gpu.py`使用Decimal60独立重算所有8784条记录。
解析附属有限核验（在课题根目录运行）：

```sh
python3 computation/audit_r016_gpu.py
python3 certificates/r016-globalization/verify_t_interface_family.py
python3 certificates/r016-globalization/audit_global_overlap_n3.py
python3 certificates/r016-globalization/verify_top2_hook_split.py
cd certificates/r016-fixed-distance
python3 fixed_d_fibre_obstruction.py
```

新位集核及P524288批次未计入本轮完成结果。

## 第26轮

行释放整数递推、独立峰集容斥、n=7 MacBook 计数、条件平均采样审计及已完成 CPU 批次，见 [本轮重放入口](../certificates/r026/REPLAY.md) 和 [报告](../rounds/026.md)。GPU 大批仍在运行，未计入本轮完成证据。

## 第27轮

n=8顶行完成数的模递推、Pratt/CRT证书、逐孔计数与已结束GPU批，见[本轮重放入口](../certificates/r027/REPLAY.md)和[报告](../rounds/027.md)。内存模型与实际RSS分别记录；未结束的后续批不计入本轮完成证据。


## 第28轮

末端折线比较、根词割集判据、n6逐孔表及n7中心四模精确计数，见[本轮重放入口](../certificates/r028/REPLAY.md)和[报告](../rounds/028.md)。本机失败模保留原因未知记录，Mac只补缺少模并核验退出。旧Mac n24批虽内部审计通过，其必要条件失败并给出>0.22的确定性二择误差界；不是可信概率证据。新CPU/GPU批留待结束后另轮核验。


## 第29轮

标签间隔终端恒等式与全层几何尾反证、n5完整尾DP、n7全部D4孔类四模CRT及已完成n24有限批，见[本轮重放入口](../certificates/r029/REPLAY.md)和[报告](../rounds/029.md)。n7矩阵的portable重放验证数学证书和收据，不冒充重新核验未归档的历史二进制。新短尾大批与GPU生产另轮处理。


## 第30轮

n7中心短尾四模Pratt/双CRT、n5完整间隔—原图几何联合DP、n3全9!插入族核验，
见[本轮重放入口](../certificates/r030/REPLAY.md)和[报告](../rounds/030.md)。
短尾原始pending-proof元数据保留，由r29已证公式解释；中心原作业已退出，新孔类
队列另行记录。正在运行的两机GPU和其余孔类不作本轮完成证据。
