# A022 第24轮归档索引

本目录是已归档的第24轮白名单材料。[`manifest-r024.json`](../../manifest-r024.json) 记录 121 个本目录
payload 及本轮入口文件的实际 SHA-256；本目录没有二进制、缓存、凭据、SSH
配置或未完成队列输出。

## 内容与证据级别

- `r024-global-normalizer/`：路径/环、叶附着恒等式、小图扫描，以及经
  审阅的度三二叉桥无限族定理。后者严格否定一般有界度图上的
  `A2 <= C(Delta) B`；不是方格原题反例。
- `r024-tree-integer/`：N=30 度三树的独立整数有限反例与实际核验。
- `r024-tree-root/`：旧 SHA `6615204c…` 的历史审计收据。匹配源码未保留，
  不可重放，不能用作冻结源的证据。
- `r024-tree-scaling/`：冻结精确树递推的已完成重跑、含祖先因子 \(G\) 的
  切边核验，以及 1440 棵 Prüfer 小树的可移植复核。
- `r024-weighted-merge/`：方格顶行臂族的严格恒等式与 n=3,4,5 整数核验；
  没有归一化远峰质量下界。
- `r024-local-distance/`、`r024-fixed-root-n16/` 和 `r024-macbook/`：已自然
  结束的有限 SMC/终端审计；它们不是渐近概率或混合认证。
- `r024-logscaled/`：最终对数缩放聚合器、夹具与回归测试；不是对 r023
  Decimal 聚合器下溢错误的修复。
- 正式中文报告见 [第24轮](../../rounds/024.md)，当前状态见 [STATE](../../STATE.md)；本目录不重复保存。

## 已执行的重放

以下均在归档前于本机实际运行并通过：

```sh
python3 -B r024-global-normalizer/verify_leaf_bridge.py
python3 -B r024-global-normalizer/verify_constants.py
python3 -B r024-tree-scaling/verify_cut_formula_archived.py
python3 -B r024-tree-scaling/verify_cut_formula.py
python3 -B r024-weighted-merge/audit_top_row_arm_mass.py
python3 -B r024-weighted-merge/verify_balanced_far_arms.py
python3 -B r024-weighted-merge/audit_top_row_branch_swaps_n3.py
python3 -B r024-logscaled/test_logscaled_sut.py
python3 -B r024-fixed-root-n16/verify_production_archive.py
(cd r024-tree-integer && python3 -B run_verify.py)
```

`r024-local-distance/ARCHIVE.md` 给出 n=16 原始记录的严格汇总重放命令。
`r024-tree-root/` 明确排除在上述重放集合外。所有仍运行的本机、GPU 或远端
队列均不在本包中。

树整数重放的计时字段随本次运行变化；不变的精确整数、源哈希及全排列小例子
与原有 `verification.stdout` 一致，新的实际输出为
`r024-tree-integer/verification-replay-r024.stdout`。

## 归档检查

归档时已按目标相对路径重建清单，并检查单文件小于 10 MiB、无二进制/凭据/SSH
配置。未纳入的运行批次不因后续完成而追溯并入本轮。
