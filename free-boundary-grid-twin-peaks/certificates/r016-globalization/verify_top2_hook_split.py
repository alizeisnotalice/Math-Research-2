#!/usr/bin/env python3
"""Exact integer audit of the right-bottom H-tree split size inequalities."""

from fractions import Fraction


def main():
    rows = {}
    for n in range(3, 129):
        N = n * n
        c = n // 2
        side = n - c - 1
        b = side * side
        a = N - b
        assert Fraction(1, 16) <= Fraction(b, N) <= Fraction(1, 4)
        assert Fraction(a * b, N * N) >= Fraction(15, 256)
        assert Fraction(2 * a * b, N * (N - 1)) >= Fraction(15, 128)
        if n in (3, 4, 5, 8, 16, 32, 64, 128):
            rows[n] = {'N': N, 'b': b, 'a': a,
                       'b_over_N': str(Fraction(b, N)),
                       'ab_over_N2': str(Fraction(a * b, N * N)),
                       'top2_interleaving_factor': str(Fraction(2 * a * b, N * (N - 1)))}
    print({'n_range': [3, 128], 'right_bottom_split_integer_bounds_verified': True,
           'samples': rows})


if __name__ == '__main__':
    main()
