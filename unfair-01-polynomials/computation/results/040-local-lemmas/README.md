# 第40轮精确局部验证

本目录保存一条条件下降定理的符号回归、两个辅助命题的精确反例，以及发现后者的有限搜索。均不是原猜想的反例。

- `check_final_low_degree_slice_closure.py`：终止系数恒等式、三次代数关系、Φ5 实分解；一般量词和完整证明见 `../../../proofs/040-low-degree-slice-descent.md`。独立审计绑定该证明 SHA-256 `9a0ad4e33eee05248deb8009a8e838878f95ac601f1e4ebef8e0cf9248855529`。
- `counterexample_exact.py`：用精确 ZZ 多项式重算 35 个不同支撑和、两个分解、gcd 及实根区间的端点符号。
- `offp_overlap_verify_phi6.py`：不导入搜索器，纯整数卷积重算 U/V/H/UV、乘积恒等式和 E13=D13=1。
- `offp_overlap_search.py` 与 `offp_overlap_p2_cap24.json`：实际有限发现范围；p=2、两侧次数≤24，Φ6 命中3次后停止。只报告实际遍历的4,052,380对，不把停止后的组合算作覆盖。

复现最终三项精确核验：

```sh
python run_checks.py
```

实际环境 Python 3.13.18 / SymPy 1.14.0；主线程从归档路径重新执行，三个子进程均返回0，输出与源码哈希见 `parent-recheck.json`。证明回归不是一般证明的替代。

发现搜索复现（另指定临时输出目录）：

```sh
python offp_overlap_search.py --p 2 --cap 24 --hit-limit 3 --output /path/to/temp/offp.json
```

论文只提供非循环共同因子的发现线索；原式定位、题录和下载哈希见 `../../../proofs/040-noncyclotomic-gcd-counterexample.md`。未公开归档论文全文。
