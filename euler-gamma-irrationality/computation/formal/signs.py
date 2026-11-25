"""Isolate positive roots of each signed kernel polynomial by Bernstein subdivision."""

import json
from math import comb, factorial, lcm
from pathlib import Path


def split(c):
    m = len(c) - 1
    left, right = [], []
    for j in range(m + 1):
        left.append(c[0] << (m - j))
        right.append(c[-1] << (m - j))
        c = [a + b for a, b in zip(c, c[1:])]
    return left, right[::-1]


def inspect(N, v):
    a = [z * lcm(*range(1, 2 * n + 1)) for n, z in enumerate(v, N)]
    while a[-1] == 0:
        a.pop()
    order = 0
    while a[0] == 0:
        order += 1
        a.pop(0)
    m = len(a) - 1
    c = [
        sum(
            a[k] * 16 ** (m - k) * comb(j, k) * (factorial(m) // comb(m, k))
            for k in range(j + 1)
        )
        for j in range(m + 1)
    ]
    stack = [(0, 0, c)]
    roots = []
    while stack:
        i, h, b = stack.pop()
        assert b[0] and b[-1] and h < 100
        signs = [x > 0 for x in b if x]
        changes = sum(x != y for x, y in zip(signs, signs[1:]))
        if changes == 1 and h >= 8:
            roots.append([i, h])
        elif changes:
            left, right = split(b)
            stack.extend([(2 * i, h + 1, left), (2 * i + 1, h + 1, right)])
    return dict(
        N=N,
        origin_order=order,
        degree=m,
        positive_at_origin=a[0] > 0,
        roots_x=sorted(roots),
    )


if __name__ == "__main__":
    d = json.loads(Path(__file__).with_name("kernel.json").read_text())[-1]
    print(json.dumps([inspect(32, r["v"]) for r in d["vectors"]], indent=2))
