# 相关性强猜想与局部交换反例

peak-association-exact.json来自原prescribed-pair DP，给出条件峰指标逐对负关联反例；root_moments_cut_ie.cpp是另一份n3..5切矩阵累计，与r019-cut-full逐项核对。中央交换见 ../../proofs/central-cut-swap-audit.md：固定较低标签峰后，有效域映射是单射，主要缺口是有效域质量。

```sh
clang++ -O3 -std=c++17 ../../computation/pair_dp.cpp -o /Volumes/data/CodexTemp/T2/a022-pair
python3 audit_peak_association.py /Volumes/data/CodexTemp/T2/a022-pair /Volumes/data/CodexTemp/T2/a022-association.json
python3 verify_cut_covariances.py
python3 audit_cut_quadratic_forms.py
python3 central_swap_witness.py
```

外切符号的均值不为零，未经中心化的E[sj(X)sj(Y)]为正不能叫正协方差。最终归档全部按中心化公式核验。最初未经中心化的解释已撤回。

复现时请先确认 data 挂载，将 TMPDIR/TMP/TEMP 及编译输出设到 `/Volumes/data/CodexTemp/T2`；Swift/Clang 模块缓存为 `/Volumes/data/CodexTemp/cache`。环境：Python 3.15，Apple Clang 21，Swift 6.4.2/Metal。源码相对路径调整不改变数学核；原运行源定位保留在本轮报告，归档哈希见 manifest。有限核验不外推到所有 n。
