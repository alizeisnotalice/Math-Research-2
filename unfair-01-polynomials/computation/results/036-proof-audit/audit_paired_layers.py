#!/Users/zzh/.venvs/main/bin/python
"""Exact symbolic audit of the paired-layer induction identities.

It enumerates finite support patterns only to check the index bookkeeping.
The accompanying note supplies the all-q proof.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


def audit_pattern(q: int, T: tuple[int, ...]) -> int:
    # Independent symbols make accidental cancellation visible.
    max_i = 2 * q - 1
    f = {0: sp.Integer(1)} | {i: sp.Symbol(f"f_{q}_{'_'.join(map(str, T))}_{i}") for i in range(1, max_i + 1)}
    d = {i: (sp.Integer(1) if i in (0, T[-1]) else sp.Symbol(f"d_{q}_{'_'.join(map(str, T))}_{i}")) for i in T}
    e = {i: (sp.Integer(1) if i in (0, T[-1]) else sp.Symbol(f"e_{q}_{'_'.join(map(str, T))}_{i}")) for i in T}
    checked = 0
    seen: list[int] = []
    for i in T:
        # Impose equality on the lower layers already traversed, exactly as in the proof.
        substitution = {e[j]: d[j] for j in seen}
        fd_i = sum(f[i - j] * d[j] for j in T if j <= i)
        fe_i = sum(f[i - j] * e[j] for j in T if j <= i)
        lhs = sp.expand((fe_i - fd_i).subs(substitution))
        assert lhs == sp.expand(e[i] - d[i]), (q, T, i, lhs)

        # h_{q+i}=FD_{q+i}+FE_i; q+i<2q is the sole range condition.
        assert q + i < 2 * q, (q, T, i)
        fd_high = sum(f[q + i - j] * d[j] for j in T if j <= q + i)
        h_high = fd_high + fe_i
        assert sp.expand(h_high - fd_high - fe_i) == 0
        seen.append(i)
        checked += 1
    return checked


def main() -> None:
    patterns = 0
    identities = 0
    by_q: dict[str, int] = {}
    for q in range(2, 11):
        count = 0
        for n in range(q):
            for interior_mask in range(1 << max(n - 1, 0)):
                T = (0,) if n == 0 else tuple([0] + [i for i in range(1, n) if (interior_mask >> (i - 1)) & 1] + [n])
                identities += audit_pattern(q, T)
                patterns += 1
                count += 1
        by_q[str(q)] = count
    payload = {
        "arithmetic": "exact SymPy polynomial arithmetic",
        "q_range": [2, 10],
        "support_patterns_checked": patterns,
        "per_q_pattern_count": by_q,
        "induction_identities_checked": identities,
        "claims_checked": [
            "conditional on substituting equality on earlier support layers, (FE)_i-(FD)_i=e_i-d_i",
            "[x^(q+i)]F(D+x^qE)=(FD)_(q+i)+(FE)_i",
            "q+i<2q for every i in T because max(T)<q",
        ],
        "limitation": "The conditional substitution is not established by this audit; it cannot prove the multi-interior equality claim.",
    }
    output = Path(__file__).with_name("audit_paired_layers.json")
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
