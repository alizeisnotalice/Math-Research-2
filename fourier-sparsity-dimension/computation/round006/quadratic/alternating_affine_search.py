#!/usr/bin/env python3
"""Exact enumeration of affine all-nonsingular pencils in Alt_4(F_2)."""
import argparse
import itertools
import json


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def rank_alt(code, n, ps):
    rows = [0] * n
    for b, (i, j) in enumerate(ps):
        if (code >> b) & 1:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    ans = 0
    for col in range(n):
        pivot = next((r for r in range(ans, n) if (rows[r] >> col) & 1), None)
        if pivot is None:
            continue
        rows[ans], rows[pivot] = rows[pivot], rows[ans]
        for r in range(n):
            if r != ans and ((rows[r] >> col) & 1):
                rows[r] ^= rows[ans]
        ans += 1
    return ans


def rref_bases(m, r):
    for piv in itertools.combinations(range(m), r):
        free = [(i, j) for i, p in enumerate(piv)
                for j in range(p + 1, m) if j not in piv]
        for bits in range(1 << len(free)):
            rows = [1 << p for p in piv]
            for b, (i, j) in enumerate(free):
                if (bits >> b) & 1:
                    rows[i] |= 1 << j
            yield tuple(rows)


def span(rows):
    a = [0]
    for x in rows:
        a += [y ^ x for y in a]
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ns = ap.parse_args()
    n = 4
    ps = pairs(n)
    m = len(ps)
    good = {x for x in range(1 << m) if rank_alt(x, n, ps) == n}
    result = {"n": n, "ambient_dimension": m, "nonsingular_count": len(good), "dimensions": {}}
    for r in range(m + 1):
        num, witness = 0, None
        for basis in rref_bases(m, r):
            sub = span(basis)
            reps = set()
            for a in range(1 << m):
                rep = min(a ^ x for x in sub)
                if rep in reps:
                    continue
                reps.add(rep)
                if all((rep ^ x) in good for x in sub):
                    num += 1
                    witness = witness or {"base": rep, "basis": list(basis), "elements": [rep ^ x for x in sub]}
        result["dimensions"][str(r)] = {"count": num, "witness": witness}
    with open(ns.out, "w") as f:
        json.dump(result, f, indent=2, sort_keys=True)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
