# H9(x³)H9(x⁵)：完整实因子认证包

目标 H9=1+x+x²+x⁸+x⁹。完整证明见 `../../../proofs/030-H935-theorem.md`。原均匀和猜想仍开放。

认证链分为：整数模分解 → 双纤维根迹/根积证明 → 循环根迹矩阵与80类完整覆盖 → 有理母根盘 → 全部80类的独立Fraction卷积。78类有严格负系数，仅恒定和标准分解保留。搜索器的96位定点负见证不被独立验证器读取。

在本目录运行（本机实际 Python 3.13.18；以下验证只依赖标准库）：

```sh
python3 verify_h9_modular.py G5-cycle-results.json
python3 verify_h9_real_layout.py
python3 verify_h935_coefficients.py --balls H9-35-root-balls.json --record /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H935-recheck.json
python3 test_h9_modular.py
python3 test_h9_semantic_tampers.py
python3 test_h9_layout_guards.py
```

`verify_h935_coefficients.py` 已包含完整Rouché前端、80唯一键覆盖、六项语义负测和三种区间角点回归。`parent-coefficient-verification.json` 是主线程从本归档包实际全量运行的完整有理上界；`parent-validation-receipt.json` 记录最终源码哈希、session和exit0。`negative-witness-summary.json` 为每类选一个负见证并向上放宽至百万分母，最大上界仍为 −683631/1000000。

根盘输入额外包含全部72个目标根的隔离。其循环部分 `(n,j)` 仅作为穷尽槽标签，数值角度未作精确绑定；系数证明不使用这些标签。循环分配由 `verify_h9_real_layout.py` 的整数指数同余独立绑定，系数由√5的严格有理界直接生成。实际系数计算仅需五个G5母盘。

生成器与搜索器（不是独立验收入口）：

```sh
python3 verify_h9_cycle_prereqs.py
python3 build_h9_real_layout.py
python3 generate_h9_roots.py
python3 search_h935.py H9-35-root-balls.json /Volumes/data/CodexTemp/T2/unfair-01-polynomials/H935-search-recheck.json
```

前两步的前提生成使用SymPy 1.14.0；根盘中心生成用mpmath，中心不论来自何种数值方法，都必须通过精确Rouché检查。完整复现时先在CodexTemp中复制本包再运行生成器，避免覆盖归档输入。生成器输出中的路径已改为相对文件名；独立整数覆盖检查重新运行通过。
