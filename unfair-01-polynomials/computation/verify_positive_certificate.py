#!/usr/bin/env python3
"""Independent exact verifier for positive-certificate-v1 JSON.

No SymPy and no import from the searcher: all polynomials are dictionaries
from exponent tuples to Fraction coefficients.  The verifier reconstructs
every convolution diagonal from I,J and separately expands both certificate
sides.  It validates only the stated finite support and certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


def clean(poly):
    return {key: value for key, value in poly.items() if value}


def add(left, right, factor=Fraction(1)):
    out = dict(left)
    for key, value in right.items():
        out[key] = out.get(key, Fraction(0))+factor*value
    return clean(out)


def scale(poly, scalar):
    return clean({key: scalar*value for key, value in poly.items()})


def mul(left, right):
    out = {}
    for a, av in left.items():
        for b, bv in right.items():
            key = tuple(x+y for x, y in zip(a, b))
            out[key] = out.get(key, Fraction(0))+av*bv
    return clean(out)


def const(n, dimension):
    return {(0,)*dimension: Fraction(n)} if n else {}


def variable(index, dimension):
    power = [0]*dimension; power[index] = 1
    return {tuple(power): Fraction(1)}


def coefficient(pair):
    return Fraction(pair[0], pair[1])


def support_check(label, support):
    if (not isinstance(support, list) or not support
            or any(type(x) is not int for x in support)
            or support[0] != 0 or any(x < 0 for x in support)
            or support != sorted(set(support))):
        raise SystemExit(
            f'{label} must be a nonempty, strictly increasing list of nonnegative integers beginning at 0')


def derived_units(I, J, diagonals):
    unit = {f'a{I[0]}', f'a{I[-1]}', f'b{J[0]}', f'b{J[-1]}'}
    for pairs in diagonals.values():
        if len(pairs) == 1:
            i, j = pairs[0]
            unit.update({f'a{i}', f'b{j}'})
    return unit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('certificate', type=Path)
    args = parser.parse_args()
    report = json.loads(args.certificate.read_text())
    if report.get('format') != 'positive-certificate-v1' or report.get('status') != 'certificate_found':
        raise SystemExit('not a successful positive-certificate-v1 file')
    I, J = report['case']['I'], report['case']['J']
    support_check('I', I)
    support_check('J', J)
    diagonals = {k: [(i, j) for i in I for j in J if i+j == k]
                 for k in sorted({i+j for i in I for j in J})}
    expected_unit = derived_units(I, J, diagonals)
    all_names = {f'a{i}' for i in I} | {f'b{j}' for j in J}
    expected_strict = all_names - expected_unit
    variables = report.get('strict_variables')
    listed_unit = report.get('unit_variables')
    if not isinstance(variables, list) or len(variables) != len(set(variables)) or set(variables) != expected_strict:
        raise SystemExit('strict_variables is not precisely the independently derived nonunit set')
    if not isinstance(listed_unit, list) or len(listed_unit) != len(set(listed_unit)) or set(listed_unit) != expected_unit:
        raise SystemExit('unit_variables is not precisely the endpoint/singleton-derived unit set')
    index = {name: pos for pos, name in enumerate(variables)}
    dimension = len(variables)
    unit = expected_unit
    def value(family, location):
        name = f'{family}{location}'
        return const(1, dimension) if name in unit else variable(index[name], dimension)
    equations = {}
    strict_witnesses = {}
    for k, pairs in diagonals.items():
        poly = {}
        for i, j in pairs:
            poly = add(poly, mul(value('a', i), value('b', j)))
        equations[k] = add(poly, const(-1, dimension))
    # Each nonunit must have its own unit contribution and a separate positive
    # support product in the same diagonal, certifying x<1 on c_k=1.
    for name in variables:
        family, location = name[0], int(name[1:])
        pairs = diagonals[location]
        own = (location, J[0]) if family == 'a' else (I[0], location)
        if own not in pairs or len(pairs) < 2:
            raise SystemExit(f'no strict own-diagonal witness for {name}')
        strict_witnesses[name] = {'diagonal': location, 'other_terms': len(pairs)-1}
    lhs = {}
    for term in report['certificate']['multipliers']:
        if term['diagonal'] not in equations:
            raise SystemExit('certificate names an absent convolution diagonal')
        poly = const(1, dimension)
        for name in term['monomial']:
            if name not in index:
                raise SystemExit('multiplier uses a non-strict or undeclared variable')
            poly = mul(poly, variable(index[name], dimension))
        lhs = add(lhs, scale(mul(poly, equations[term['diagonal']]), coefficient(term['coefficient'])))
    rhs = {}
    atom_positive = []
    for atom in report['certificate']['positive_atoms']:
        positive, one_minus = atom['positive'], atom['one_minus']
        if (len(positive) != len(set(positive)) or len(one_minus) != len(set(one_minus))
                or set(positive) & set(one_minus)
                or not set(positive).issubset(index) or not set(one_minus).issubset(index)):
            raise SystemExit('positive atom is not a squarefree product over declared strict variables')
        poly = const(1, dimension)
        for name in positive:
            poly = mul(poly, variable(index[name], dimension))
        for name in one_minus:
            poly = mul(poly, add(const(1, dimension), variable(index[name], dimension), Fraction(-1)))
        q = coefficient(atom['coefficient'])
        if q < 0:
            raise SystemExit('negative positive-atom coefficient')
        if q:
            atom_positive.append(atom)
        rhs = add(rhs, scale(poly, q))
    if not atom_positive or lhs != rhs:
        raise SystemExit('exact polynomial equality failed')
    print(json.dumps({
        'status': 'verified', 'variables': variables,
        'equation_count': len(equations), 'multiplier_terms': len(report['certificate']['multipliers']),
        'positive_atoms': len(atom_positive), 'strict_witnesses': strict_witnesses,
        'claim_boundary': 'The certificate is strictly positive only on the declared open cube 0<x<1; it does not classify arbitrary supports or search bases.'
    }, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
