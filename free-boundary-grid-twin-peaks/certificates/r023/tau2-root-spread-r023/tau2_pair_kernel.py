#!/usr/bin/env python3
"""Exact top-two completion kernels h({x,y}) on P_n square, n=3,4,5.

Uses the independent-set inclusion-exclusion identity from top-two-inclusion-
exclusion.md.  All arithmetic is Python integers; no sampling is used.
"""
from collections import defaultdict
from itertools import combinations
import json
from math import factorial
from fractions import Fraction


def independent_sets(n):
    rows = [r for r in range(1 << n) if not (r & (r << 1))]
    out = []
    def rec(y, prev, mask):
        if y == n:
            out.append(mask)
            return
        for row in rows:
            if not (row & prev):
                rec(y + 1, row, mask | (row << (y * n)))
    rec(0, 0, 0)
    return out


def audit(n):
    N = n * n
    full = (1 << N) - 1
    adj = [0] * N
    for v in range(N):
        x, y = v % n, v // n
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            xx, yy = x + dx, y + dy
            if 0 <= xx < n and 0 <= yy < n:
                adj[v] |= 1 << (yy * n + xx)
    def neigh(mask):
        ans = 0
        q = mask
        while q:
            bit = q & -q
            q -= bit
            ans |= adj[bit.bit_length() - 1]
        return ans
    sets = independent_sets(n)
    sets.sort(key=int.bit_count)
    C = {0: factorial(N)}
    for s in sets[1:]:
        total = 0
        q = s
        while q:
            bit = q & -q
            q -= bit
            total += C[s ^ bit]
        den = s.bit_count() + neigh(s).bit_count()
        assert total % den == 0
        C[s] = total // den
    assert len(C) == len(sets)
    a_root = [0] * N
    for s, value in C.items():
        if not s:
            continue
        sign = 1 if s.bit_count() & 1 else -1
        q = s
        while q:
            bit = q & -q
            q -= bit
            a_root[bit.bit_length() - 1] += sign * value
    assert all(z > 0 for z in a_root)

    h = {}
    for x, y in combinations(range(N), 2):
        if (adj[x] >> y) & 1:
            continue
        closed = (1 << x) | (1 << y) | adj[x] | adj[y]
        total = 0
        for s in sets:
            if not (s & closed):
                total += -C[s] if s.bit_count() & 1 else C[s]
        assert total % (N * (N - 1)) == 0
        h[x, y] = total // (N * (N - 1))

    # Fixed roots representative of corner, noncorner boundary, and centre.
    roots = [(0, 0), (0, 1), (n // 2, n // 2)]
    root_rows = {}
    fixed_root_distance_failures = []
    half_tail = []
    for xy in roots:
        x = xy[1] * n + xy[0]
        vals = []
        for y in range(N):
            if y == x or (adj[x] >> y) & 1:
                continue
            e = (min(x, y), max(x, y))
            d = abs(x % n - y % n) + abs(x // n - y // n)
            vals.append({'y': [y % n, y // n], 'distance': d, 'h': h[e]})
        vals.sort(key=lambda z: (z['distance'], z['y']))
        root_rows[str(xy)] = vals
        for u in vals:
            for v in vals:
                if u['distance'] < v['distance'] and u['h'] > v['h']:
                    fixed_root_distance_failures.append({'x': list(xy), 'near': u, 'far': v})
    # All fixed roots, for the concrete half-side-scale diagnostic d>n/2.
    for x in range(N):
        values = [(abs(x % n - y % n) + abs(x // n - y // n), h[(min(x, y), max(x, y))])
                  for y in range(N) if y != x and not ((adj[x] >> y) & 1)]
        den = sum(v for _, v in values)
        num = sum(v for d, v in values if 2 * d > n)
        half_tail.append({'x': [x % n, x // n], 'numerator': num, 'denominator': den,
                          'fraction': str(Fraction(num, den)),
                          'a_x': a_root[x], 'H_x_over_a_x': str(Fraction(den, a_root[x]))})
    # h is D4 invariant; explicit scalar summaries simplify inspection.
    by_dist = defaultdict(list)
    for (x, y), value in h.items():
        d = abs(x % n - y % n) + abs(x // n - y // n)
        by_dist[d].append(value)
    return {
        'n': n, 'independent_sets': len(sets), 'nonadjacent_pairs': len(h),
        'A1': sum(a_root),
        'h_distance_ranges': {str(d): {'minimum': min(v), 'maximum': max(v),
                                        'sum': sum(v), 'pairs': len(v)}
                              for d, v in sorted(by_dist.items())},
        'representative_fixed_root_rows': root_rows,
        'fixed_root_distance_failure': (fixed_root_distance_failures[0]
                                        if fixed_root_distance_failures else None),
        'strict_half_side_tail_by_fixed_high_root': half_tail,
        'minimum_strict_half_side_tail': min(half_tail, key=lambda z: Fraction(z['numerator'], z['denominator'])),
        'maximum_H_over_a': max(half_tail, key=lambda z: Fraction(z['denominator'], z['a_x'])),
    }


print(json.dumps({'status': 'passed', 'audits': [audit(n) for n in (3, 4, 5)]},
                 indent=2, sort_keys=True))
