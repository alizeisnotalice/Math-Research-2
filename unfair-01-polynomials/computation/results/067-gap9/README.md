# 第 9 行内部卷积的精确排除

[完整定理与证明](../../../proofs/067-ninth-interior-row.md)。

全 27 形中有 5 形直接用第 8 行上界排除。其余 22 形逐项用结式、不可约分解、Sturm 根计数核验，两个本原因子通过整性与实根筛后都被有理区间负号排除。另一个词典序 Gröbner 入口对这 22 形产生相同本原因子。使用 Python 3.13.18、SymPy 1.14.0。

    python -B verify_gap9.py --out /Volumes/data/CodexTemp/T2/a30-gap9-exact.json
    python -B probe_gap9.py --out /Volumes/data/CodexTemp/T2/a30-gap9-groebner.jsonl

独立脚本将同一卷积模型直接构造为多项式，重算有理 Sturm 根计数并用未约简分式核负号。

    python -B independent/audit_gap9.py verification.json --out /Volumes/data/CodexTemp/T2/a30-gap9-independent.json
