# 混合替换正性下降的精确交叉核验

[完整解析定理](../../../proofs/032-crossed-descent-degree-bound.md)：对互素整数1<p<q，在H与H(x^p)H(x^q)二元、FG=H且两个交叉因子非负的完整前提下，degH≤4q时F,G必非负。一般定理由解析论证承担，以下有限测试核查公式和反例模型。

主线程实际运行：

```sh
python3 verify_descent_identities.py
python3 verify_binary_H_nonbinary_C.py
/Users/zzh/.venvs/main/bin/python verify_cross_recurrences.py
```

第一项以独立稀疏整数卷积核验10对素数、两种显式非负恒等式、余数格以及等号4q边界的两项/三项系数式。第二项以另一套稠密整数卷积验证36对素数，特别是H仍二元而C_pq=2的精确例。第三项用SymPy形式系数验证两端递推及p3q5、p3q7两个等号分支。所有成功运行退出0；`parent-validation.json` 与 `symbolic-run.json` 记录实际环境、源码及输出哈希。

环境澄清：本机默认 `python3` 实为3.15.0，仅上述标准库程序可直接运行；它未安装SymPy，首次符号检查导入失败。随后使用已有Python3.13.18数学环境（SymPy1.14.0）完成检查，没有安装新库。第31轮文字统一写3.13.18，而其主线程摘要环境字段记录3.15.0；本轮另明确用3.13.18全量复跑H10的17类Fraction验证，输出与第31轮逐字节一致。本轮不将环境摘要混作每项历史运行的解释器证据。

公式分工和实际核验范围见 [残数网格证明](../../../proofs/032-residue-grid-obstruction.md)。此计算包不提供任意次数的P31证明，也未解决一般均匀和猜想。
