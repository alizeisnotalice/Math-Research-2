# 解析引理的有限核验

对应 `../../proofs/` 中 interface-surface-cost-obstruction、cut-peak-budget、detachable-branch-insertion、detachable-branch-root-image、weighted-root-covariance。

在本目录运行（输出重定向到临时卷）：

```sh
python3 audit_surface_ratios.py
python3 audit_peak_budget.py
python3 audit_detachable_branch_n3.py
python3 audit_first_merge_root_n4.py 3
python3 audit_first_merge_root_n4.py 4
clang++ -O3 -std=c++17 ../../computation/pair_dp.cpp -o /Volumes/data/CodexTemp/T2/a022-r018-pair-dp
python3 audit_pair_reflection.py /Volumes/data/CodexTemp/T2/a022-r018-pair-dp /Volumes/data/CodexTemp/T2/a022-r018-pair.json
```

整数枚举与解析证明分开保存。浮点展示只用于显示整数比。部分脚本写当前目录同名小 JSON，复现时可复制脚本至临时卷运行。所有矩形统一接口常数界被细长梯图反例否定；不能据此外推到有界纵横比。

本轮只研究自由边界方格、原图距离、均匀标号条件 K=2。有限数据不证明渐近结论。

复现环境为 Python 3.15、Apple Clang 21、Swift 6.4.2 / Metal。所有编译与输出请放 `/Volumes/data/CodexTemp/T2`，Clang/Swift 模块缓存使用 `/Volumes/data/CodexTemp/cache`。
