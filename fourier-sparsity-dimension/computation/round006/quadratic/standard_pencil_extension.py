#!/usr/bin/env python3
"""Exact one-direction extension test for the standard quadratic pencil, k=3.

The pencil has polar matrices [[A, I+N],[I+N^T,0]] with A alternating and
N strictly upper triangular.  We test every remaining alternating direction D
for whether adjoining D preserves nonsingularity on all 64 pencil points.
"""
import argparse
import json


def pair_index(n):
    ps = [(i, j) for i in range(n) for j in range(i + 1, n)]
    return ps, {p: b for b, p in enumerate(ps)}


def add_entry(code, lookup, i, j):
    if i == j:
        return code
    if i > j:
        i, j = j, i
    return code ^ (1 << lookup[(i, j)])


def rank_alt(code, n, ps):
    rows = [0] * n
    for b, (i, j) in enumerate(ps):
        if (code >> b) & 1:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    r = 0
    for c in range(n):
        p = next((i for i in range(r, n) if (rows[i] >> c) & 1), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(n):
            if i != r and ((rows[i] >> c) & 1):
                rows[i] ^= rows[r]
        r += 1
    return r


def span(basis):
    ans = [0]
    for x in basis:
        ans += [y ^ x for y in ans]
    return ans


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    k, n = 3, 6
    ps, lookup = pair_index(n)
    base = 0
    basis = []
    for i in range(k):
        base = add_entry(base, lookup, i, k + i)
    # A-block alternating directions.
    for i in range(k):
        for j in range(i + 1, k):
            basis.append(add_entry(0, lookup, i, j))
    # Strict-upper cross-block directions N_{ij}, i<j.
    for i in range(k):
        for j in range(i + 1, k):
            basis.append(add_entry(0, lookup, i, k + j))
    direction = set(span(basis))
    pencil = {base ^ x for x in direction}
    good = {x for x in range(1 << len(ps)) if rank_alt(x, n, ps) == n}
    extend = []
    for d in range(1 << len(ps)):
        if d in direction:
            continue
        if all((x ^ d) in good for x in pencil):
            extend.append(d)
    result = {
        "k": k, "ambient_alternating_dimension": len(ps),
        "pencil_direction_dimension": len(basis), "pencil_size": len(pencil),
        "all_pencil_nonsingular": pencil <= good,
        "one_direction_extensions": extend,
        "extension_count": len(extend),
    }
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2, sort_keys=True)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
