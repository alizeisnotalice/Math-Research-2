#!/usr/bin/env python3
"""Search a bounded low-degree strict-positive convolution certificate.

This is a discovery tool.  It fixes endpoint and singleton-diagonal variables
to one, forms e_k=c_k-1, and solves a rational-recovered LP in the ansatz
  sum p_{k,m} m e_k = sum q_{P,N} prod(P) prod(1-N), q>=0.
The basis bounds are only a search rule; failure is reported as unknown.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import linprog
import sympy as S


CASE = {
    'name': 'I0247_J0246891012131518',
    'I': [0, 2, 4, 7],
    'J': [0, 2, 4, 6, 8, 9, 10, 12, 13, 15, 18],
}


def diagonals(I, J):
    return {k: [(i, j) for i in I for j in J if i+j == k]
            for k in sorted({i+j for i in I for j in J})}


def monomials(symbols, degree):
    out = [S.Integer(1)]
    for r in range(1, degree+1):
        out.extend(S.prod(parts) for parts in combinations_with_replacement(symbols, r))
    return out


def polynomial_vector(poly, symbols, basis):
    P = S.Poly(S.expand(poly), *symbols)
    return [P.coeff_monomial(m) for m in basis]


def encode_fraction(value):
    q = Fraction(value)
    return [q.numerator, q.denominator]


def variable_name(family, index):
    return f'{family}{index}'


def validate_supports(I, J):
    """Require the normalized, finite supports used by this searcher."""
    for label, support in (('I', I), ('J', J)):
        if (not support or any(type(x) is not int for x in support)
                or support[0] != 0 or any(x < 0 for x in support)
                or support != sorted(set(support))):
            raise ValueError(
                f'{label} must be a nonempty, strictly increasing list of nonnegative integers beginning at 0')


def derived_units(I, J, D=None):
    """Variables forced to one by endpoint or singleton convolution equations."""
    D = diagonals(I, J) if D is None else D
    units = {('a', I[0]), ('a', I[-1]), ('b', J[0]), ('b', J[-1])}
    for pairs in D.values():
        if len(pairs) == 1:
            i, j = pairs[0]
            units.update({('a', i), ('b', j)})
    return units


def encode_monomial(m):
    """JSON's repeated-name encoding is unambiguous also for powers a2**r."""
    if m == 1:
        return []
    powers = m.as_powers_dict()
    encoded = []
    for symbol in sorted(powers, key=str):
        exponent = powers[symbol]
        if not exponent.is_Integer or exponent < 0:
            raise ValueError(f'unsupported multiplier monomial {m!r}')
        encoded.extend([str(symbol)]*int(exponent))
    return encoded


def discover(I, J, multiplier_degree, atom_degree, max_denominator):
    validate_supports(I, J)
    D = diagonals(I, J)
    units = derived_units(I, J, D)
    keys = [('a', x) for x in I] + [('b', x) for x in J]
    strict_keys = [key for key in keys if key not in units]
    symbol_result = S.symbols(' '.join(variable_name(*key) for key in strict_keys)) if strict_keys else ()
    symbols = ((symbol_result,) if len(strict_keys) == 1 else tuple(symbol_result))
    value = {key: S.Integer(1) for key in units}
    value.update({key: sym for key, sym in zip(strict_keys, symbols)})
    equations = {k: S.expand(sum(value[('a', i)]*value[('b', j)] for i, j in pairs)-1)
                 for k, pairs in D.items()}
    equations = {k: e for k, e in equations.items() if e != 0}
    # With no strict variables the model is already completely decided.
    # A collision is an immediate contradiction; otherwise all coefficients
    # are one and there is no nonuniform point to eliminate.
    if not strict_keys:
        collisions = {k: len(pairs) for k, pairs in D.items() if len(pairs) > 1}
        status = 'unit_collision_infeasible' if collisions else 'trivially_uniform_no_strict_variables'
        return {
            'format': 'positive-certificate-v1', 'case': {'name': CASE['name'], 'I': I, 'J': J},
            'strict_variables': [], 'unit_variables': sorted(variable_name(*key) for key in units),
            'basis': {'multiplier_total_degree_at_most': multiplier_degree,
                      'positive_atom_total_degree_at_most': atom_degree,
                      'atom_form': 'product(positive) * product(1-minus)'},
            'matrix': {'polynomial_rows': 0, 'normalization_rows': 0,
                       'free_multiplier_columns': 0, 'positive_atom_columns': 0,
                       'lp_nonnegative_columns': 0, 'elapsed_seconds': 0.0},
            'solver': {'name': 'not run (all variables forced to one)', 'success': True,
                       'message': status},
            'status': status,
            'unit_collision_diagonals': collisions,
        }
    multipliers = monomials(symbols, multiplier_degree)
    product_columns = [(k, m, S.expand(m*e)) for k, e in equations.items() for m in multipliers]
    atoms = []
    for r in range(1, atom_degree+1):
        for picks in combinations(range(len(symbols)), r):
            for signs in product((0, 1), repeat=r):
                positive = tuple(str(symbols[t]) for t, sign in zip(picks, signs) if sign == 0)
                one_minus = tuple(str(symbols[t]) for t, sign in zip(picks, signs) if sign == 1)
                atom = S.prod(symbols[t] if sign == 0 else 1-symbols[t]
                              for t, sign in zip(picks, signs))
                atoms.append((positive, one_minus, S.expand(atom)))
    max_degree = max(S.Poly(p, *symbols).total_degree()
                     for _, _, p in product_columns + [(None, None, p) for _, _, p in atoms])
    basis = monomials(symbols, max_degree)
    L = S.Matrix.hstack(*[S.Matrix(polynomial_vector(p, symbols, basis)) for _, _, p in product_columns])
    R = S.Matrix.hstack(*[S.Matrix(polynomial_vector(p, symbols, basis)) for _, _, p in atoms])
    n_left, n_atom = L.shape[1], R.shape[1]
    matrix = np.array(S.Matrix.hstack(L, -L, -R), dtype=float)
    Aeq = np.vstack([matrix, np.r_[np.zeros(2*n_left), np.ones(n_atom)]])
    beq = np.r_[np.zeros(matrix.shape[0]), 1.0]
    started = perf_counter()
    # Minimize multiplier L1 norm after normalizing the nonnegative RHS cone.
    answer = linprog(np.r_[np.ones(2*n_left), np.zeros(n_atom)], A_eq=Aeq, b_eq=beq,
                     bounds=(0, None), method='highs')
    elapsed = perf_counter()-started
    report = {
        'format': 'positive-certificate-v1', 'case': {'name': CASE['name'], 'I': I, 'J': J},
        'strict_variables': [str(s) for s in symbols],
        'unit_variables': sorted(variable_name(*key) for key in units),
        'basis': {'multiplier_total_degree_at_most': multiplier_degree,
                  'positive_atom_total_degree_at_most': atom_degree,
                  'atom_form': 'product(positive) * product(1-minus)'},
        'matrix': {'polynomial_rows': matrix.shape[0], 'normalization_rows': 1,
                   'free_multiplier_columns': n_left, 'positive_atom_columns': n_atom,
                   'lp_nonnegative_columns': 2*n_left+n_atom, 'elapsed_seconds': elapsed},
        'solver': {'name': 'SciPy HiGHS linprog', 'success': bool(answer.success), 'message': answer.message},
    }
    if not answer.success:
        report['status'] = 'unknown_within_declared_basis'
        return report
    rational = [S.Rational(Fraction(float(z)).limit_denominator(max_denominator)) for z in answer.x]
    p = [rational[q]-rational[n_left+q] for q in range(n_left)]
    q = rational[2*n_left:]
    lhs = S.expand(sum(coef*poly for coef, (_, _, poly) in zip(p, product_columns)))
    rhs = S.expand(sum(coef*poly for coef, (_, _, poly) in zip(q, atoms)))
    exact = S.expand(lhs-rhs) == 0 and all(z >= 0 for z in q) and any(z > 0 for z in q)
    report['certificate'] = {
        'multipliers': [
            {'diagonal': k,
             'monomial': encode_monomial(m),
             'coefficient': encode_fraction(coef)}
            for coef, (k, m, _) in zip(p, product_columns) if coef != 0
        ],
        'positive_atoms': [
            {'positive': list(pos), 'one_minus': list(neg), 'coefficient': encode_fraction(coef)}
            for coef, (pos, neg, _) in zip(q, atoms) if coef != 0
        ],
        'exact_recovery_check': bool(exact),
    }
    report['status'] = 'certificate_found' if exact else 'floating_candidate_not_exactly_recovered'
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--multiplier-degree', type=int, default=1)
    parser.add_argument('--atom-degree', type=int, default=2)
    parser.add_argument('--I', help='comma-separated support, default is the built-in case')
    parser.add_argument('--J', help='comma-separated support, default is the built-in case')
    parser.add_argument('--name', help='case name for the certificate JSON')
    parser.add_argument('--max-denominator', type=int, default=4096)
    args = parser.parse_args()
    if args.multiplier_degree < 0 or args.atom_degree < 1 or args.max_denominator < 1:
        raise SystemExit('invalid degree bound')
    if bool(args.I) != bool(args.J):
        raise SystemExit('--I and --J must be supplied together')
    I = [int(x) for x in args.I.split(',')] if args.I else CASE['I']
    J = [int(x) for x in args.J.split(',')] if args.J else CASE['J']
    try:
        validate_supports(I, J)
    except ValueError as exc:
        raise SystemExit(str(exc))
    report = discover(I, J, args.multiplier_degree, args.atom_degree, args.max_denominator)
    report['case']['name'] = args.name or report['case']['name']
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'status': report['status'], 'matrix': report['matrix'], 'solver': report['solver'],
                      'certificate_terms': len(report.get('certificate', {}).get('multipliers', [])),
                      'positive_atoms': len(report.get('certificate', {}).get('positive_atoms', []))}, ensure_ascii=False))


if __name__ == '__main__':
    main()
