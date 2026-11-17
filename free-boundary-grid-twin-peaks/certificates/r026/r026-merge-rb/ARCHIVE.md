# 首次合并 Rao--Blackwell 审计

`audit_merge_rao_blackwell.py` 以 Fraction 穷举 3×3 的全部 (9!) 标号，并对 32128
个固定首次合并 fibre 逐一验证：对根交换不变、交错不变的 (H=1,d^2,F_0,G_0)，
完整双峰和等于 (\tau=2\) 的 (L H\) 和 (mH\) 和，其中
\(m=q(q-1)/(2ab)\)。它特别记录 67200 个 (L\ne m\) 的词，故该结论不是逐词替换。

最小 r026 GPU trace（n=3、P=17）由 `verify_merge_rb_trace.py` 独立重放。精确 fibre
审计默认从同一 `certificates/` 根下既有 r025 `r025-rank-bias/` 读取其最小旧 trace，
以避免重复归档；必要时用 `A022_CERTIFICATES_ROOT` 指定该根。两种 trace 都仅核验有限
实现路径，不能认证 SMC 混合、方差改善、normalizer 精度或原问题。二进制、较大 trace、
生产批均未收录。
