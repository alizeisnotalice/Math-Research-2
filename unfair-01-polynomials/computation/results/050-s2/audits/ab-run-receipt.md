# S2 AB 全批运行回执

本回执记录本 Agent 的实际 AB 树运行；没有重跑父方的屏幕或 HGQ 验证器。

## 冻结输入

来源目录：

```text
/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round050-parent/s2-screen-complete
```

`snapshot_s2_screen.py` 实际检查 `core-000.json` 至 `core-358.json` 均为
`status=complete`、索引匹配，逐源记录原始 SHA-256、`binding` 和幸存行数，拒绝
重复八元组。所得冻结输入：

```text
inputs_s2_all359_frozen.json
SHA-256 4f85bb7b53217381f7d4e49cc8d9f8efd204a059021490912dc0974a2b6b1661
359 source cores; 9,775 distinct rows
```

父屏幕的完整摘要报告原始枚举 37,808,085、Newton 后 1,477,286、正轴后
845,726、外 Schur 后 13,212、内 Schur 幸存 9,775；这些是必要母式筛的范围，
不是本 AB 运行自身重新声明的结论。

## 实际 AB 命令和结果

在本 Agent 的 zsh 中实际运行：

```sh
nice -n 10 python3 \
  /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round050-s2-bfs/s2_exact_prefix_bfs.py \
  /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round050-s2-bfs/inputs_s2_all359_frozen.json \
  --out /Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round050-s2-bfs/bfs_s2_all359_depth96.json \
  --depth 96 --state-cap 25000 --seconds-per-core 8 --quiet
```

解释器为 `/Users/zzh/.local/bin/python3`（Python 3.15.0，Clang 23.1.3）。
输出摘要：9,775 条 `dead_all_prefixes`，无 `unknown_*`，最晚死亡层 67，最大活态
数 4。输出和程序哈希分别为：

```text
bfs_s2_all359_depth96.json
b3ae9b3091dd739012c43077f401740da30ab4092de12c870b95bb0346fb2fa2

s2_exact_prefix_bfs.py
eb418a74a7d0e230f872b4798042361ed43d51d0f8e5c949ef11db3b5f9cd90a
```

树在每一行精确计算 \(\mathbb Z[a]/(a^2-Ta+N)\)。它允许任意符号、任意有限
次数的整数 (G,Q)，只加入已证必要条件 (G_3=Q_3=0)、(A_5=B_5=0)。
任意有限余因子最高次后补零，因而对应母式树空时可严格排除它；未把 96 当作
余因子次数上限。

## 独立 HGQ 比对

`compare_s2_ab_hgq.py` 实际读取父目录
`round050-parent/s2-hgq-complete/core-*.json` 的全部 359 文件。它逐 core 检查
源 SHA 与母式集合，并逐行要求 status、完整逐层 trace（含 live、state SHA 和
三类 reject 计数）完全相同。

```text
s2_ab_hgq_comparison_all359.json
SHA-256 2ac93c6d7d67752f5901da1fb135b99ce426461da0f8a13e4bd6afea32392dbb
359 cores; 9,775 exact traces; 0 discrepancy
```

该紧凑账本为每母式保存 source core/源 SHA、AB 输入与结果 SHA、AB 记录索引、
死亡层、最大活态数和双方 trace SHA。完整 AB trace 输出仅保留在临时目录，可由
冻结输入重建。

## 回归和退出

先前实际执行的 smoke 回归包括有理二元控制、一个人工非有理短前缀，以及全部
395 个既有低端前缀逐态对照（359 个不同核心），均通过。完整 AB 命令已退出；
本回执不保留后台进程。
