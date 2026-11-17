#!/usr/bin/env python3
"""Exact finite audits; no claim that enumeration proves all dimensions."""
import argparse
import hashlib
import itertools
import json
import platform
import time
from pathlib import Path


def basis(vectors):
    pivots = {}
    for x in vectors:
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                break
    return tuple(pivots.values())


def span(vectors):
    s = {0}
    for x in vectors:
        s |= {v ^ x for v in tuple(s)}
    return frozenset(s)


def matrix_rows(mask, n):
    rows = [0] * n
    for k, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        if (mask >> k) & 1:
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
    return rows


def wedge(a, b, n):
    return sum(((((a >> i) & 1) * ((b >> j) & 1)) ^
                 (((a >> j) & 1) * ((b >> i) & 1))) << k
               for k, (i, j) in enumerate(itertools.combinations(range(n), 2)))


def check_family(points, planes, n, counts):
    a = min(points)
    differences = [x ^ a for x in points]
    d = len(points).bit_length() - 1
    assert len(points) == 1 << d
    assert len(basis(differences)) == d
    common = set(planes[a]) - {0}
    for x in points:
        common.intersection_update(planes[x])
    core_dim = len(basis(v for x in points for v in planes[x]))
    rank4 = any(len(basis(matrix_rows(x, n))) == 4 for x in differences)
    if common:
        typ = "star"
    else:
        typ = "bounded_core"
        assert core_dim <= 4 and d <= 3, (n, sorted(points), core_dim, d)
    if rank4:
        assert core_dim == 4 and d <= 3 and not common
    key = f"d{d}_{typ}_rank4diff{int(rank4)}"
    counts[key] = counts.get(key, 0) + 1


def all_subspaces(s):
    levels = [{frozenset({0})}]
    while levels[-1]:
        nxt = set()
        for u in levels[-1]:
            for x in s - u:
                coset = frozenset(v ^ x for v in u)
                if coset <= s:
                    nxt.add(u | coset)
        if not nxt:
            break
        levels.append(nxt)
    return levels


def audit(n, full_cosets):
    ncoords = n * (n - 1) // 2
    ranktwo = {x for x in range(1 << ncoords)
               if len(basis(matrix_rows(x, n))) == 2}
    planes = {x: span(matrix_rows(x, n)) for x in ranktwo}
    # Independent construction of every decomposable nonzero alternating form.
    decomposable = {wedge(a, b, n) for a in range(1, 1 << n)
                    for b in range(1, 1 << n) if a != b}
    assert decomposable == ranktwo
    counts = {}
    if full_cosets:
        levels = all_subspaces(frozenset(range(1 << ncoords)))
        affine_tested = 0
        for level in levels:
            for u in level:
                remaining = set(range(1 << ncoords))
                while remaining:
                    a = min(remaining)
                    points = frozenset(a ^ v for v in u)
                    remaining.difference_update(points)
                    affine_tested += 1
                    if points <= ranktwo:
                        check_family(points, planes, n, counts)
        scope = "all affine subspaces of the entire alternating-form space"
    else:
        anchor = min(ranktwo)
        allowed = frozenset(x ^ anchor for x in ranktwo)
        levels = all_subspaces(allowed)
        affine_tested = sum(map(len, levels))
        for level in levels:
            for u in level:
                check_family(frozenset(x ^ anchor for x in u), planes, n, counts)
        scope = "all constant-rank-two affine families containing the fixed form e0 wedge e1"
    return {"n": n, "rank_two_forms": len(ranktwo), "scope": scope,
            "subspace_level_counts": [len(x) for x in levels],
            "affine_tested": affine_tested, "valid_classifications": counts}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=4)
    p.add_argument("--all-cosets", action="store_true")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    t = time.monotonic()
    result = audit(a.n, a.all_cosets)
    result.update({"seconds": time.monotonic() - t, "python": platform.python_version(),
                   "platform": platform.platform(), "arithmetic": "exact GF(2), integer masks",
                   "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   "failed_assertions": 0})
    Path(a.out).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
