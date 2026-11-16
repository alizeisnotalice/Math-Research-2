#!/usr/bin/env python3
"""Integer audit of corrected internal-level bookkeeping for the constant 60."""
from fractions import Fraction


def main():
    levels = 0
    for q in range(1, 4097):
        J = (q - 1).bit_length() - 1
        if q == 1:
            assert J == -1
            continue
        assert 2 ** J < q <= 2 ** (J + 1)
        for ell in range(J + 1):
            cell = (q + 2 ** ell - 1) // (2 ** ell)
            assert cell * 2 ** ell <= 2 * q
            levels += 1
        s0 = sum((Fraction(1, 2 ** t) for t in range(J + 1)), Fraction())
        s1 = sum((Fraction(t, 2 ** t) for t in range(J + 1)), Fraction())
        assert s0 <= 2 and s1 <= 2
    print({'q_range': [1, 4096], 'corrected_internal_levels_checked': levels,
           'q_equals_one_has_no_internal_node': True,
           'dyadic_geometric_coefficients_at_most_two': True})


if __name__ == '__main__':
    main()
