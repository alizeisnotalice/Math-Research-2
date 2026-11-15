#!/usr/bin/env python3
"""Find a rational Farkas dual for one fixed positive-certificate template.

This intentionally lives outside the frozen discovery/verifier pair.  If it
returns a certificate, its functional annihilates every allowed multiplier
column and is at least one on every allowed positive atom.  Hence no nonzero
nonnegative combination of those atoms can equal an allowed left side.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
import sympy as S

from positive_certificate_search import (  # noqa: E402
    derived_units, diagonals, encode_monomial, monomials, polynomial_vector,
    validate_supports, variable_name,
)


def atoms_for(symbols, atom_degree):
    atoms = []
    for r in range(1, atom_degree+1):
        for picks in combinations(range(len(symbols)), r):
            for signs in product((0, 1), repeat=r):
                positive = tuple(str(symbols[t]) for t, sign in zip(picks, signs) if sign == 0)
                one_minus = tuple(str(symbols[t]) for t, sign in zip(picks, signs) if sign == 1)
                atom = S.prod(symbols[t] if sign == 0 else 1-symbols[t]
                              for t, sign in zip(picks, signs))
                atoms.append((positive, one_minus, S.expand(atom)))
    return atoms


def full_quadratic_atoms(symbols):
    """All one/two literal products, including repeated-variable factors."""
    literals = [(str(s), False, s) for s in symbols] + [(str(s), True, 1-s) for s in symbols]
    atoms = []
    for degree in (1, 2):
        for chosen in combinations_with_replacement(range(len(literals)), degree):
            positive = tuple(literals[t][0] for t in chosen if not literals[t][1])
            one_minus = tuple(literals[t][0] for t in chosen if literals[t][1])
            atom = S.prod(literals[t][2] for t in chosen)
            atoms.append((positive, one_minus, S.expand(atom)))
    return atoms


def fraction_pair(x):
    return [x.numerator, x.denominator]


def search(I, J, multiplier_degree, atom_degree, max_denominator, full_quadratic=False):
    validate_supports(I, J)
    D = diagonals(I, J)
    units = derived_units(I, J, D)
    keys = [('a', x) for x in I] + [('b', x) for x in J]
    strict = [key for key in keys if key not in units]
    if not strict:
        raise ValueError('this dual template needs at least one strict variable')
    symbol_result = S.symbols(' '.join(variable_name(*key) for key in strict))
    symbols = (symbol_result,) if len(strict) == 1 else tuple(symbol_result)
    value = {key: S.Integer(1) for key in units}
    value.update({key: sym for key, sym in zip(strict, symbols)})
    equations = {
        k: S.expand(sum(value[('a', i)]*value[('b', j)] for i, j in pairs)-1)
        for k, pairs in D.items()
    }
    equations = {k: e for k, e in equations.items() if e != 0}
    multipliers = monomials(symbols, multiplier_degree)
    columns = [(k, m, S.expand(m*e)) for k, e in equations.items() for m in multipliers]
    atoms = full_quadratic_atoms(symbols) if full_quadratic else atoms_for(symbols, atom_degree)
    max_degree = max(S.Poly(p, *symbols).total_degree()
                     for _, _, p in columns + [(None, None, p) for _, _, p in atoms])
    basis = monomials(symbols, max_degree)
    L = S.Matrix.hstack(*[S.Matrix(polynomial_vector(p, symbols, basis)) for _, _, p in columns])
    R = S.Matrix.hstack(*[S.Matrix(polynomial_vector(p, symbols, basis)) for _, _, p in atoms])
    # y=y_plus-y_minus, y^T L=0, y^T atom >= 1 for every permitted atom.
    l_float, r_float = np.array(L, dtype=float), np.array(R, dtype=float)
    rows = L.shape[0]
    answer = linprog(
        np.ones(2*rows),
        A_ub=-np.hstack([r_float.T, -r_float.T]), b_ub=-np.ones(R.shape[1]),
        A_eq=np.hstack([l_float.T, -l_float.T]), b_eq=np.zeros(L.shape[1]),
        bounds=(0, None), method='highs')
    report = {
        'format': ('positive-template-farkas-dual-fullquadratic-v1' if full_quadratic
                   else 'positive-template-farkas-dual-v1'),
        'case': {'I': I, 'J': J},
        'strict_variables': [str(s) for s in symbols],
        'unit_variables': sorted(variable_name(*key) for key in units),
        'basis': {'multiplier_total_degree_at_most': multiplier_degree,
                  'positive_atom_total_degree_at_most': atom_degree,
                  'atom_form': ('all products with replacement from x and (1-x)'
                                if full_quadratic else 'product(positive) * product(1-minus)')},
        'matrix': {'polynomial_rows': rows, 'left_columns': L.shape[1],
                   'atom_columns': R.shape[1], 'lp_nonnegative_columns': 2*rows},
        'solver': {'name': 'SciPy HiGHS linprog', 'success': bool(answer.success),
                   'message': answer.message},
    }
    if not answer.success:
        report['status'] = 'dual_not_found'
        return report
    y = [S.Rational(Fraction(float(answer.x[t]-answer.x[rows+t])).limit_denominator(max_denominator))
         for t in range(rows)]
    exact_left = all(v == 0 for v in (L.T*S.Matrix(y)))
    atom_values = list(R.T*S.Matrix(y))
    exact_atoms = all(v >= 1 for v in atom_values)
    report['functional'] = [
        {'monomial': encode_monomial(m), 'coefficient': fraction_pair(Fraction(c))}
        for m, c in zip(basis, y) if c != 0
    ]
    report['exact_check'] = {
        'annihilates_all_multiplier_columns': bool(exact_left),
        'minimum_atom_value': fraction_pair(Fraction(min(atom_values))) if atom_values else None,
        'all_atom_values_at_least_one': bool(exact_atoms),
    }
    report['status'] = 'farkas_dual_found' if exact_left and exact_atoms else 'floating_dual_not_exactly_recovered'
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--I', required=True)
    parser.add_argument('--J', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--multiplier-degree', type=int, default=1)
    parser.add_argument('--atom-degree', type=int, default=2)
    parser.add_argument('--max-denominator', type=int, default=1000000)
    parser.add_argument('--full-quadratic', action='store_true',
                        help='allow x^2, x(1-x), (1-x)^2 as degree-two atoms')
    args = parser.parse_args()
    I = [int(x) for x in args.I.split(',')]
    J = [int(x) for x in args.J.split(',')]
    if args.full_quadratic and args.atom_degree != 2:
        raise SystemExit('--full-quadratic currently has exactly degree-two atoms')
    report = search(I, J, args.multiplier_degree, args.atom_degree, args.max_denominator,
                    full_quadratic=args.full_quadratic)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': report['status'], 'matrix': report['matrix'],
                      'solver': report['solver'], 'exact_check': report.get('exact_check')}))


if __name__ == '__main__':
    main()
