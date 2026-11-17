# 已作废的首次条件账本

`core6_linked_conditional_ledger.json` 是 2026-11-17 的首次运行产物，**不得使用**。

原因是脚本当时把

\[
\left\lceil\frac{-24-b^3}{3}\right\rceil
\]

的 Python 整数下取整写错，导致 `ell` 的下端错误；该错误不满足预期的
\((b,\ell)\mapsto(-b,-\ell)\) 对称性，且实际漏掉了正 `b` 的可行行。

修复后的脚本将输出写到
`core6_linked_conditional_ledger_corrected.json`。只有该新文件可用于后续工作，且它仍然仅是已审核根界条件下的有限账本，并非排除证明。
