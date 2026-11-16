#!/usr/bin/env python3
"""Small exact audits of the T-interface matching construction."""

from itertools import combinations, permutations
from math import comb, factorial
from random import Random


def neighbors(v, n):
    x, y = v
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        q = (x + dx, y + dy)
        if 0 <= q[0] < n and 0 <= q[1] < n:
            yield q


def births(word, n):
    seen, count = set(), 0
    for v in word:
        if not any(q in seen for q in neighbors(v, n)):
            count += 1
        seen.add(v)
    return count


def regions(m):
    n = 2 * m + 1
    A = {(x, y) for x in range(m) for y in range(n - 1)}
    B = {(x, y) for x in range(m + 1, n) for y in range(n - 1)}
    C = {(m, y) for y in range(n)} | {(x, n - 1) for x in range(n)}
    w, c = (m, n - 2), (m, n - 1)
    Cp = C - {w, c}
    parent = {(m, y): (m - 1, y) for y in range(n - 2)}
    parent |= {(x, n - 1): (x, n - 2) for x in range(m)}
    parent |= {(x, n - 1): (x, n - 2) for x in range(m + 1, n)}
    assert set(parent) == Cp
    assert len(set(parent.values())) == len(Cp)
    return n, A, B, Cp, parent, w, c


def extensions(alpha, leaves):
    k, M = len(leaves), len(alpha)
    gates = alpha[M - k:]
    tail = gates + leaves
    good = []
    for candidate in permutations(tail):
        at = {v: i for i, v in enumerate(candidate)}
        if (all(at[gates[j]] < at[gates[j + 1]] for j in range(k - 1))
                and all(at[gate] < at[leaf] for gate, leaf in zip(gates, leaves))):
            good.append(candidate)
    return good


def audit_n3():
    n, A, B, Cp, parent, w, c = regions(1)
    a_orders, b_orders = list(permutations(sorted(A))), list(permutations(sorted(B)))
    words = set()
    for ao in a_orders:
        for bo in b_orders:
            for a_slots in combinations(range(4), 2):
                alpha, ia, ib = [], 0, 0
                for i in range(4):
                    if i in a_slots:
                        alpha.append(ao[ia]); ia += 1
                    else:
                        alpha.append(bo[ib]); ib += 1
                if {alpha[0], alpha[1]} != {ao[0], bo[0]}:
                    continue
                leaves = [v for _, v in sorted((alpha.index(parent[v]), v) for v in Cp)]
                for tail in extensions(alpha, leaves):
                    for remainder in ((w, c), (c, w)):
                        word = tuple(alpha[:1] + list(tail) + list(remainder))
                        assert births(word, n) == 2
                        words.add(word)
    e3 = factorial(6) // (2**3 * factorial(3))
    expected = 2 * len(a_orders) * len(b_orders) * (2 * comb(2, 1)) * e3
    assert len(words) == expected
    return {"n": n, "words": len(words), "E3": e3}


def audit_n5():
    n, A, B, Cp, parent, w, c = regions(2)
    M, k = len(A) + len(B), len(Cp)
    assert (M, k) == (16, 7)
    e7 = factorial(14) // (2**7 * factorial(7))
    assert e7 == 135135
    rng = Random(20261116)
    for _ in range(200):
        alpha = list(A | B)
        rng.shuffle(alpha)
        rank = {v: i + 1 for i, v in enumerate(alpha)}
        leaves = [v for _, v in sorted((rank[parent[v]], v) for v in Cp)]
        for j, leaf in enumerate(leaves, 1):
            assert rank[parent[leaf]] <= M - k + j
    revealed = A | B | Cp
    assert all(any(q in revealed for q in neighbors(v, n)) for v in (w, c))
    return {"n": n, "E7": e7, "gate_samples": 200, "tail_orders": 2}


if __name__ == "__main__":
    print(audit_n3())
    print(audit_n5())
