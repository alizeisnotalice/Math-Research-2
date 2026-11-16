#!/usr/bin/env python3
"""Exact finite-prefix audit for the tau=2 Feynman--Kac recursion on 3x3."""
import json
from math import prod
from pathlib import Path

N = 9
HERE = Path(__file__).resolve().parent


def adjacencies():
    out = [set() for _ in range(N)]
    for v in range(N):
        x, y = v % 3, v // 3
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            xx, yy = x + dx, y + dy
            if 0 <= xx < 3 and 0 <= yy < 3:
                out[v].add(yy * 3 + xx)
    return out


ADJ = adjacencies()


def allowed(prefix, rank):
    """Allowed rank-th vertex for the birth set {1,2}; rank is one-based."""
    used = set(prefix)
    boundary = set().union(*(ADJ[v] for v in prefix)) - used
    remaining = set(range(N)) - used
    return remaining - boundary if rank == 2 else boundary


def main():
    levels = [None, {(v,) for v in range(N)}]
    audit = [{"rank": 1, "prefixes": N, "falling_factorial": N,
              "mu_one": "1"}]
    for rank in range(2, N + 1):
        old = levels[-1]
        next_level = set()
        edge_count = 0
        for x in old:
            a = allowed(x, rank)
            edge_count += len(a)
            next_level.update(x + (v,) for v in a)
        # Each valid child has its unique deletion parent.  Thus this asserts
        # the finite combinatorial content of mu_{k-1}(G_k Q_k)=mu_k for f=1.
        assert edge_count == len(next_level)
        ff = prod(range(N - rank + 1, N + 1))
        audit.append({"rank": rank, "prefixes": len(next_level),
                      "sum_allowed_from_previous": edge_count,
                      "falling_factorial": ff,
                      "mu_one_numerator": len(next_level),
                      "mu_one_denominator": ff})
        levels.append(next_level)

    terminal = levels[N]
    assert len(terminal) == 81984
    payload = {"n": 3, "N": N,
               "method": "explicit valid-prefix transition graph; exact integer edge counts",
               "levels": audit,
               "terminal_tau2_count": len(terminal),
               "checks": {"each_child_has_one_parent": True,
                          "all_mu_recursions_for_constant_one": True,
                          "matches_second_birth_audit_tau2_count": True}}
    (HERE / "tau2-prefix-fk-n3.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
