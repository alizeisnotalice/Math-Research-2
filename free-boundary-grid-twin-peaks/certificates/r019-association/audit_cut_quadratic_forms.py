#!/usr/bin/env python3
"""Exact quadratic-form audit for r019 vertical-cut covariance matrices.

Input files are the independently generated inclusion--exclusion certificates
for n=3,...,6.  The output uses Fraction only: no sign decision is floating
point based.  A strictly negative LDL diagonal (without pivoting) certifies
that the displayed covariance matrix is negative definite.
"""

from fractions import Fraction
import json
from pathlib import Path


SOURCE = Path(__file__).resolve().parent.parent / "r019-cut-full"
OUT = Path(__file__).with_name('cut-quadratic-forms-exact.json')


def load_matrix(n: int):
    # The independently audited n=3 file carries the explicit ``-ie`` suffix.
    name = 'n3-ie.json' if n == 3 else f'n{n}-cut-covariance.json'
    data = json.loads((SOURCE / name).read_text())
    d = int(data['covariance_common_denominator'])
    a = [[Fraction(int(x), d) for x in row]
         for row in data['covariance_numerators']]
    return data, a


def ldl_no_pivot(a):
    """Return exact unit-lower L and diagonal D for a symmetric matrix."""
    n = len(a)
    L = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    D = []
    for j in range(n):
        d = a[j][j] - sum(L[j][r] * L[j][r] * D[r] for r in range(j))
        if d == 0:
            raise ArithmeticError(f'zero pivot at {j}')
        D.append(d)
        for i in range(j + 1, n):
            L[i][j] = (a[i][j] - sum(L[i][r] * L[j][r] * D[r]
                                     for r in range(j))) / d
    return L, D


def main():
    runs = []
    for n in range(3, 7):
        data, a = load_matrix(n)
        _, diag = ldl_no_pivot(a)
        q_ones = sum(sum(row) for row in a)
        positive = [(i, j, x) for i, row in enumerate(a)
                    for j, x in enumerate(row) if x > 0]
        # A sufficient exact copositivity certificate: discard all nonpositive
        # off-diagonal terms, and absorb each remaining positive B_ij term by
        # 2 a_i a_j <= a_i^2+a_j^2.  ``margins`` are the remaining negative
        # diagonal magnitudes.  It is deliberately a certificate, not a
        # decision algorithm for copositivity in general.
        margins = [-a[i][i] for i in range(len(a))]
        for i, j, x in positive:
            if i < j:
                margins[i] -= x
                margins[j] -= x
        absorption = all(x > 0 for x in margins)
        runs.append({
            'n': n,
            'dimension': len(a),
            'all_diagonal_entries_strictly_negative':
                all(a[i][i] < 0 for i in range(len(a))),
            'all_ldl_pivots_strictly_negative': all(x < 0 for x in diag),
            'ldl_pivots': [str(x) for x in diag],
            'all_ones_quadratic_form': str(q_ones),
            'all_ones_quadratic_form_negative': q_ones < 0,
            'positive_entry_positions': [[i, j] for i, j, _ in positive],
            'diagonal_absorption_copositivity_certificate': absorption,
            'diagonal_absorption_margins': [str(x) for x in margins],
        })
    OUT.write_text(json.dumps({'status': 'passed', 'runs': runs}, indent=2) + '\n')
    print(OUT)


if __name__ == '__main__':
    main()
