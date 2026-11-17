# 六次 linked `{2,3}` 核：条件账本

本目录按父线程给出的、**尚待 section 独立审核**的根环/牛顿界，精确枚举

\[
P=1+b x+a x^2+(ba+\ell)x^3+d x^4+e x^5+x^6,
\qquad a^2-Ta+N=0,
\]

的母参数 `(b,T,N,ell)`。其中 `b` 非零整数，`a` 是实非有理二次代数整数；`d,e` 尚未展开。

脚本只用整数平方比较同时核验两嵌入的绝对值约束，避免浮点根。它报告各筛选阶段的完整计数，并把剩余行写到 `core6_linked_conditional_ledger_corrected.json`。`d,e` 的理论盒大小为 `43*15=645`，因此输出还报告若以后展开该盒的精确候选量。首次错误输出已在 `INVALIDATED_FIRST_RUN.md` 中明确作废。

**状态：条件性预处理，不是排除证明，也不是完整的有限枚举证书。** 在根界经独立审核前，JSON 中的行仅为父线程给定不等式下的精确结果。

复现：

```sh
/Users/zzh/.venvs/main/bin/python generate_core6_linked_ledger.py
```
