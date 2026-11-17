# 第39轮精确切片核验

主线程从本目录实际运行五个脚本，回执为 `parent-rerun.json`。均为 Python 3.13.18；代数程序使用 SymPy 1.14.0，互素例程序只用标准库整数算术。

```sh
python check_coprime_slice_descent.py
python check_fixed_point_without_S.py
python check_counterexample.py
python check_low_degree_S_closure.py
python check_common_factor_c_boundary.py
```

分别验证双长侧 primitive 公平例的二元乘积和完整 F₂ Euclid 互素性、指定57未知系统的精确满列秩、以及不动点过强推论的实代数边界反例。前两者不是任意模型的自动性质；第三者不满足完整交叉假设，不是原猜想反例。一般定理另见 proofs/039-coprime-slices-and-rigidity.md。
