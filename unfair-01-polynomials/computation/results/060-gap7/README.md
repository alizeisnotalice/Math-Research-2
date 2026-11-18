# 第 7 项缺口的次数下界

在 \(A_1,B_1>0,\ C_7=0\) 下，任一因子的次数至少为 20。另一因子的次数和系数域的扩张次数任意。

- [次数 0 至 15](../../../proofs/060-c7-minimum-factor-degree.md)
- [次数 16](../../../proofs/060-c7-degree16-mod163.md)
- [次数 17、18](../../../proofs/060-c7-degrees-17-18.md)
- [次数 19](../../../proofs/060-c7-degree-19.md)

在本目录使用 Python 3.13.18 重放。第一条命令依赖 SymPy 1.14.0，其他三个入口用 Python 标准库。以包外未存在的路径替换新输出文件名。

    python3 -B certify_c7_small_degrees.py --out /Volumes/data/CodexTemp/T2/c7-small-degrees-new.json
    python3 -B verify_gap7_degree16.py --out /Volumes/data/CodexTemp/T2/c7-degree16-new.json
    python3 -B verify_degrees_17_18.py --out /Volumes/data/CodexTemp/T2/c7-degrees-17-18-new.json
    python3 -B verify_degree_19.py --out /Volumes/data/CodexTemp/T2/c7-degree19-new.json

回执保存在同目录。次数 17、18 和 19 两套脚本对 640、512 个比特组合分别精确重算；次数 19 的计数包含提前拒绝的右右分支的代数余式检查。证明使用的是施加实数正性以后的子集。
