#!/usr/bin/env python3
"""Enumerate necessary short-G support patterns that evade mod-p isolation.

For deg(G)<2q, the proved q-class pairing and g_0,g_q>0 give
Gsupport = T union (T+q), where 0 is in T subset {0,...,q-1}.
Since this support contains q-separated pairs, binary C forces no p-separated
pair in Gsupport.  Equivalently T has neither difference p nor q-p.

Any occupied singleton residue class of Gsupport modulo p isolates every f_i
in A=F(x^p)G(x^q), hence rules out a signed F.  We enumerate only patterns
with no singleton; they are necessary support-level survivors, not examples.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import gcd
from pathlib import Path


def enumerate_T(p: int, q: int):
    """Yield all T containing 0 with no differences p or q-p."""
    forbidden_diffs = {p, q - p}
    chosen = {0}

    def dfs(x: int):
        if x == q:
            yield tuple(sorted(chosen))
            return
        if all(x - d not in chosen for d in forbidden_diffs):
            chosen.add(x)
            yield from dfs(x + 1)
            chosen.remove(x)
        yield from dfs(x + 1)

    yield from dfs(1)


def support(T: tuple[int, ...], q: int) -> tuple[int, ...]:
    return tuple(sorted(set(T) | {x + q for x in T}))


def injective_check(S: tuple[int, ...], p: int, q: int) -> bool:
    sums: set[int] = set()
    for a in S:
        for b in S:
            v = p * a + q * b
            if v in sums:
                return False
            sums.add(v)
    return True


def non_singleton(S: tuple[int, ...], p: int) -> bool:
    counts = [0] * p
    for x in S:
        counts[x % p] += 1
    return all(c != 1 for c in counts)


def one_pair(p: int, q: int) -> dict:
    total = 0
    survivors: list[dict] = []
    for T in enumerate_T(p, q):
        total += 1
        S = support(T, q)
        # Independent direct check of the support-addition criterion.
        assert injective_check(S, p, q)
        if non_singleton(S, p):
            survivors.append({"T": list(T), "G_support": list(S)})
    return {
        "p": p,
        "q": q,
        "T_patterns_binary_C_compatible": total,
        "no_mod_p_singleton_survivors": survivors,
        "no_mod_p_singleton_count": len(survivors),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-q", type=int, default=20)
    ap.add_argument("--q-ge-2p-only", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    for q in range(3, args.max_q + 1):
        for p in range(2, q):
            if gcd(p, q) != 1 or (args.q_ge_2p_only and q < 2 * p):
                continue
            rows.append(one_pair(p, q))
    result = {
        "schema": "short-G-singleton-support-v1",
        "range": {"q<=": args.max_q, "coprime_1<p<q": True, "q_ge_2p_only": args.q_ge_2p_only},
        "scope": "Necessary support patterns only. A survivor is not a real factorization or a counterexample.",
        "summary": {
            "pairs": len(rows),
            "binary_C_compatible_T_patterns": sum(x["T_patterns_binary_C_compatible"] for x in rows),
            "no_singleton_survivors": sum(x["no_mod_p_singleton_count"] for x in rows),
            "pairs_with_survivor": sum(bool(x["no_mod_p_singleton_survivors"]) for x in rows),
        },
        "rows": rows,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.out.write_text(encoded)
    print(json.dumps(result["summary"], sort_keys=True))
    print("output_sha256=" + hashlib.sha256(encoded.encode()).hexdigest())
    print("source_sha256=" + hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
