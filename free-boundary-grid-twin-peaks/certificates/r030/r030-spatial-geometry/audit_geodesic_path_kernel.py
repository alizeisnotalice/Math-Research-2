#!/usr/bin/env python3
"""Exact finite audit for the geodesic induced-path root kernel.

Small paths are checked by enumerating all label orders.  The L-shaped
geodesic path uses exact binomial/Fraction arithmetic for its barycentre and
variance; it illustrates why endpoint midpoints are not interchangeable with
the binomial barycentre.
"""

from fractions import Fraction
from itertools import permutations
from math import comb


def adjacent_path(i, j):
    return abs(i - j) == 1


def one_peak_roots_by_permutations(vertices):
    count = [0] * len(vertices)
    total = 0
    for word in permutations(vertices):  # descending labels
        seen = {word[0]}
        valid = True
        for v in word[1:]:
            if not any(adjacent_path(v, u) for u in seen):
                valid = False
                break
            seen.add(v)
        if valid:
            count[word[0]] += 1
            total += 1
    return total, count


def binomial_path_moments(coords):
    L = len(coords) - 1
    weights = [comb(L, i) for i in range(L + 1)]
    den = 1 << L
    mux = sum(Fraction(weights[i] * coords[i][0], den) for i in range(L + 1))
    muy = sum(Fraction(weights[i] * coords[i][1], den) for i in range(L + 1))
    var = sum(
        Fraction(weights[i], den) * ((coords[i][0] - mux) ** 2 + (coords[i][1] - muy) ** 2)
        for i in range(L + 1)
    )
    return (mux, muy), var


def l_path(m):
    """(0,0)->(m,0)->(m,m), an induced geodesic path of length 2m."""
    return [(i, 0) for i in range(m + 1)] + [(m, j) for j in range(1, m + 1)]


def main():
    print('small full-permutation checks')
    for L in range(1, 7):
        total, roots = one_peak_roots_by_permutations(tuple(range(L + 1)))
        target = [comb(L, i) for i in range(L + 1)]
        assert roots == target
        assert total == 1 << L
        print(f'  path_vertices={L + 1}: A1={total}, roots={roots}')

    print('L-shaped induced geodesic paths')
    for m in (2, 3, 4, 8, 16):
        coords = l_path(m)
        (mx, my), var = binomial_path_moments(coords)
        alpha = sum(Fraction(comb(2 * m, i) * (i - m), 1 << (2 * m))
                    for i in range(m + 1, 2 * m + 1))
        assert (mx, my) == (m - alpha, alpha)
        assert var <= Fraction(2 * m, 4)
        # E[(I-m)_+]^2 <= E[(I-m)^2] = m/2, recorded without floating point.
        assert alpha * alpha <= Fraction(m, 2)
        print(f'  m={m}: alpha={alpha}, mu=({mx},{my}), var={var}, L/4={Fraction(m,2)}')


if __name__ == '__main__':
    main()
