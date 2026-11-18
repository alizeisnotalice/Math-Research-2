# 第47轮复现说明

下列命令从本目录执行，临时输出指定在 data 卷的空目录内。实测 Python 3.15.0、Apple clang 21.0.0。完备证明、有限精确计数和有限采样的适用范围见第47轮报告。

## 路径修复与原降序词

```sh
python3 -B band-gap/verify.py | cmp - band-gap/verification.json
python3 -B band-gap/audit_path_moment.py | cmp - band-gap/audit_path_moment.stdout
python3 -B erase-rank-gap/audit_rank_gap.py | cmp - erase-rank-gap/audit_rank_gap.stdout
```

比较输出包含实际 Python 版本时，其他版本可逐项对照其余字段。脚本中的种子、整数或有理数构造与全部尺寸是此次实际门槛。

## 九阶账本和保护峰分母

```sh
python3 -B n9-family/verify.py "$scratch/family.json"
cmp "$scratch/family.json" n9-family/summary.json
python3 -B n9-upper/combine_family_bound.py n9-upper/summary.json \
  n9-family/summary.json "$scratch/conditional.json"
cmp "$scratch/conditional.json" n9-upper/original-n9-lower.json
```

家族账本保存三个原批次的摘要、全部90项原模结果文本及其回执。首批6项来自第46轮，两个新批各42项。重放复核全部输入哈希、场状态层数、模数和整数恢复。若重新计算模结果，按 `../r046/layered/REPLAY.md` 构建原计数器，依 `residue-ledger.json` 三个批次各自的 `planned` 选择区域。`local-area40` 采用3工作进程，每进程512 MiB的表预算，合计1536 MiB。`macbook-area41` 采用3工作进程，1024 MiB每项，合计3072 MiB；前两机器角色由该账本标明。表预算与实际最大RSS分别记录。

分母的独立小规模校准使用完整 \(9!\) 枚举和六个小区域的全排列。原生产还读已知 \(A_2(5)\)、\(A_2(7)\) 核验新上界。公开脚本把该两项已知计数作为固定输入 `n9-upper/calibration-counts.json`，其余计数逻辑保持一致。将该目录的 `count_upper.py`、`protected_mod.cpp`、`pratt-61bit.json`、`calibration-counts.json` 和 `reference_crt.py` 拷到空的 scratch，再运行

```sh
clang++ -std=c++20 -O3 -Wall -Wextra -Werror protected_mod.cpp -o protected_mod
python3 -B count_upper.py actual --workers 12
```

资源回执随机器不同。应比较 `upper_bounds`、`regions`、`counter_sha256`、`n3_full_permutations_checked` 与 `primes`，本轮公开重编译的60项计数均与原生产逐项相同。

## Metal 有限诊断

```sh
python3 -B n20-metal/verify.py n20-metal/n20.jsonl n20-metal/receipt.json \
  n20-metal/resource.txt \
  ../r026/r026-merge-rb-calibration/source/tau2_merge_rb_calibration_r026.swift \
  | cmp - n20-metal/summary.json
```

完整CLI参数、种子域、源码与原始结果哈希及退出记录均在账本中。从所指 Swift 源码按第26轮方式编译后，可用 `receipt.json` 的参数重跑。重新运行的计时依设备而变。权重稳定性诊断不提供混合证明或置信保证。
