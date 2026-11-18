# 最大平方因子的局部检验

[任意次数定理、边界惯性及提升证明](../../../proofs/065-square-factor-ramification.md) · [独立审查](independent/audit-zh.md)。

验证使用 Python 3.13.18、SymPy 1.14.0。平方自由七次参数的 5,760 个合成算术检查和重复参数的两种簇测试分别保存为 `square-boundary.json`、`exceptional-lift.json`。合成多项式用于检验模运算。

```sh
python -B verify_square_boundary.py --out /Volumes/data/CodexTemp/T2/a30-square-boundary-new.json
python -B verify_exceptional_lift.py --out /Volumes/data/CodexTemp/T2/a30-exceptional-lift-new.json
python -B check_c7_product.py c7-fixture.json
```

通用检查入口接收升幂排列的二元系数 JSON 列表，要求给定缺口分支的前十二项。每处返回平方因子次数，并在恰为 7 时重建 (J) 和各不可约余因子的提升等式。带重根时要求相关值为素域非零平方。通过必要条件的目标报告为未决。小输入 `c7-fixture.json` 在两素数下均由平方次数为零排除。
