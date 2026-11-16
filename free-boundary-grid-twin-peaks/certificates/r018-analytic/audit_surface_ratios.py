#!/usr/bin/env python3
"""Exact small-rectangle diagnostics for corner-root/interface ratios.

This is an exact integer subset-DP.  It is deliberately limited to at most
16 vertices in the default run: it is an audit, not an asymptotic estimator.
"""

from math import comb, factorial, log


def grid(w, h):
    V = [(x, y) for y in range(h) for x in range(w)]
    pos = {v: i for i, v in enumerate(V)}
    nbr = []
    for x, y in V:
        mask = 0
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            if (x + dx, y + dy) in pos:
                mask |= 1 << pos[x + dx, y + dy]
        nbr.append(mask)
    return V, nbr


def B(nbr, roots):
    N = len(nbr); full = (1 << N) - 1
    rootmask = sum(1 << x for x in roots)
    dp = [0] * (1 << N); dp[0] = 1
    for S in range(1 << N):
        z = dp[S]
        if not z:
            continue
        rest = full ^ S
        while rest:
            bit = rest & -rest; rest ^= bit
            v = bit.bit_length() - 1
            if bit & rootmask or nbr[v] & S:
                dp[S | bit] += z
    return dp[full]


def audit(w, h):
    V, nbr = grid(w, h); N = w * h
    roots = [B(nbr, (i,)) for i in range(N)]
    A1 = sum(roots)
    I = tuple(i for i, (x, y) in enumerate(V) if x == w - 1)
    BI = B(nbr, I)
    # Disjoint decomposition by the sole possible peak outside I.
    T = BI
    for q in range(N):
        if q not in I:
            T += B(nbr, I + (q,)) - BI
    assert T >= BI >= 1
    corner = roots[0]
    # The two rigorous ladder inequalities proved in REPORT.md are checked
    # here only as a guard against a convention mismatch in this DP.
    if w == 2:
        assert T >= factorial(h) ** 2
        assert A1 <= h * factorial(h) * (2 ** (2 * h - 1)) * comb(2 * h, h)
    return {
        'shape': f'{w}x{h}', 'A1': A1, 'corner_root_A1': corner,
        'max_root_A1': max(roots),
        'max_over_corner': f'{max(roots)}/{corner}',
        'corner_fraction': f'{corner}/{A1}',
        'minus_log_corner_fraction': -log(corner / A1),
        'B_right': BI, 'B_right_over_A1': f'{BI}/{A1}',
        'T_right': T, 'T_over_A1': f'{T}/{A1}',
        'log_T_over_A1': log(T / A1),
        'perimeter': 2 * (w + h),
        'all_counts_exact': True,
    }


if __name__ == '__main__':
    for w, h in ((2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7),
                 (2, 8), (3, 3), (3, 4), (4, 4)):
        print(audit(w, h))
