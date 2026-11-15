#!/usr/bin/env python3
"""Exhaustive n=3 audit of A2(beta>0)/A1 = E_1 T.

All 9! label orders are inspected.  The left side uses the deletion formula
M>b directly, while the right side uses the proposed one-peak order statistic
T=sum_{a=1}^{N-2} a 1{v_{a+2} has no neighbour in S_a}.  Thus neither side
calls the other's first-merge implementation.
"""
import itertools
import json
from pathlib import Path

n = 3
N = n * n
ADJ = [set() for _ in range(N)]
for v in range(N):
    x, y = v % n, v // n
    for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
        if 0 <= X < n and 0 <= Y < n:
            ADJ[v].add(Y * n + X)


def peaks(order):
    pos = [0] * N
    for i, v in enumerate(order):
        pos[v] = i  # lower i = higher label
    return [v for v in range(N) if all(pos[v] < pos[x] for x in ADJ[v])]


def beta_positive(order, ps):
    # p[0] is the global maximum; p[1] is the lower peak.
    u, v = ps
    pos = [0] * N
    for i, z in enumerate(order):
        pos[z] = i
    M = min(pos[x] for x in ADJ[v])  # rank of highest v-neighbour
    restricted = [z for z in range(N) if z != v]
    new_peaks = [z for z in restricted if all(x == v or pos[z] < pos[x] for x in ADJ[z])]
    b = min((pos[z] for z in new_peaks if z != u), default=N)  # rank N represents label 0
    # Labels decrease as order ranks increase, so M>b in label notation iff rank(M)<rank(b).
    return M < b


def T(order):
    total = 0
    for a in range(1, N - 1):  # a=1,...,N-2
        prior = set(order[:a])
        v = order[a + 1]
        if not (ADJ[v] & prior):
            total += a
    return total


def main():
    a1 = support = sum_T = 0
    for order in itertools.permutations(range(N)):
        ps = peaks(order)
        if len(ps) == 1:
            a1 += 1
            sum_T += T(order)
        elif len(ps) == 2:
            ps.sort(key=order.index)
            support += beta_positive(order, ps)
    assert a1 == 29568
    assert support == 79360
    assert sum_T == support
    out = {"n": n, "A1": a1, "beta_positive_count": support,
           "sum_onepeak_T": sum_T, "identity_verified": True}
    target = Path(__file__).parent / "results/support-T-n3-r006.json"
    target.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
