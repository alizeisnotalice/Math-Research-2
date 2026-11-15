#!/usr/bin/env python3
"""Independent Fraction verifier for positive-template-farkas-dual-v1.

It deliberately does not import the search code.  It rebuilds supports,
endpoint/singleton units, all affine-multiplied convolution columns, and the
complete squarefree positive-atom family from the JSON metadata.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path


def clean(poly):
    return {key: value for key, value in poly.items() if value}


def add(left, right, factor=Fraction(1)):
    out = dict(left)
    for key, value in right.items():
        out[key] = out.get(key, Fraction(0)) + factor*value
    return clean(out)


def mul(left, right):
    out = {}
    for a, av in left.items():
        for b, bv in right.items():
            key = tuple(x+y for x, y in zip(a, b))
            out[key] = out.get(key, Fraction(0)) + av*bv
    return clean(out)


def const(value, dimension):
    return {(0,)*dimension: Fraction(value)} if value else {}


def variable(index, dimension):
    exponent = [0]*dimension
    exponent[index] = 1
    return {tuple(exponent): Fraction(1)}


def pair(value):
    if not isinstance(value, list) or len(value) != 2 or type(value[0]) is not int or type(value[1]) is not int:
        raise SystemExit('coefficient must be an integer numerator/denominator pair')
    if value[1] == 0:
        raise SystemExit('zero denominator')
    return Fraction(value[0], value[1])


def support_check(label, support):
    if (not isinstance(support, list) or not support
            or any(type(x) is not int for x in support)
            or support[0] != 0 or any(x < 0 for x in support)
            or support != sorted(set(support))):
        raise SystemExit(f'{label} must be nonempty, increasing, nonnegative and start at zero')


def diagonals(I, J):
    return {k: [(i, j) for i in I for j in J if i+j == k]
            for k in sorted({i+j for i in I for j in J})}


def derived_units(I, J, diagonal_map):
    result = {f'a{I[0]}', f'a{I[-1]}', f'b{J[0]}', f'b{J[-1]}'}
    for pairs in diagonal_map.values():
        if len(pairs) == 1:
            i, j = pairs[0]
            result.update({f'a{i}', f'b{j}'})
    return result


def exponent_for(names, index, dimension):
    exponent = [0]*dimension
    for name in names:
        if name not in index:
            raise SystemExit('functional uses a non-strict or undeclared variable')
        exponent[index[name]] += 1
    return tuple(exponent)


def evaluate(functional, polynomial):
    return sum(functional.get(exponent, Fraction(0))*coefficient
               for exponent, coefficient in polynomial.items())


def multiplier_monomials(dimension, degree):
    out = [()]
    for size in range(1, degree+1):
        out.extend(choice for choice in combinations_with_replacement(range(dimension), size))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dual', type=Path)
    args = parser.parse_args()
    report = json.loads(args.dual.read_text())
    formats = {
        'positive-template-farkas-dual-v1': 'squarefree',
        'positive-template-farkas-dual-fullquadratic-v1': 'fullquadratic',
    }
    atom_mode = formats.get(report.get('format'))
    if atom_mode is None or report.get('status') != 'farkas_dual_found':
        raise SystemExit('not a successful supported Farkas-dual file')
    I, J = report['case']['I'], report['case']['J']
    support_check('I', I); support_check('J', J)
    diagonal_map = diagonals(I, J)
    expected_unit = derived_units(I, J, diagonal_map)
    expected_strict = ({f'a{i}' for i in I} | {f'b{j}' for j in J}) - expected_unit
    variables, units = report.get('strict_variables'), report.get('unit_variables')
    if not isinstance(variables, list) or len(variables) != len(set(variables)) or set(variables) != expected_strict:
        raise SystemExit('strict-variable partition differs from independently reconstructed partition')
    if not isinstance(units, list) or len(units) != len(set(units)) or set(units) != expected_unit:
        raise SystemExit('unit-variable partition differs from independently reconstructed partition')
    basis_spec = report.get('basis', {})
    multiplier_degree = basis_spec.get('multiplier_total_degree_at_most')
    atom_degree = basis_spec.get('positive_atom_total_degree_at_most')
    if type(multiplier_degree) is not int or multiplier_degree < 0 or type(atom_degree) is not int or atom_degree < 1:
        raise SystemExit('invalid template degree bounds')
    expected_atom_form = ('all products with replacement from x and (1-x)'
                          if atom_mode == 'fullquadratic' else 'product(positive) * product(1-minus)')
    if basis_spec.get('atom_form') != expected_atom_form:
        raise SystemExit('unexpected atom family declaration')
    if atom_mode == 'fullquadratic' and atom_degree != 2:
        raise SystemExit('full-quadratic verifier supports exactly degree-two atoms')
    index, dimension = {name: pos for pos, name in enumerate(variables)}, len(variables)
    maximum_degree = max(multiplier_degree+2, atom_degree)
    functional = {}
    for item in report.get('functional', []):
        exponent = exponent_for(item.get('monomial', []), index, dimension)
        if sum(exponent) > maximum_degree or exponent in functional:
            raise SystemExit('invalid or duplicated functional monomial')
        functional[exponent] = pair(item.get('coefficient'))
    if not functional:
        raise SystemExit('empty dual functional')
    def value(family, location):
        name = f'{family}{location}'
        return const(1, dimension) if name in expected_unit else variable(index[name], dimension)
    equations = {}
    for k, pairs in diagonal_map.items():
        poly = {}
        for i, j in pairs:
            poly = add(poly, mul(value('a', i), value('b', j)))
        equations[k] = add(poly, const(-1, dimension))
    annihilated = 0
    for equation in equations.values():
        for term in multiplier_monomials(dimension, multiplier_degree):
            monomial = const(1, dimension)
            for coordinate in term:
                monomial = mul(monomial, variable(coordinate, dimension))
            if evaluate(functional, mul(monomial, equation)) != 0:
                raise SystemExit('functional does not annihilate an allowed multiplier column')
            annihilated += 1
    atom_values = []
    if atom_mode == 'squarefree':
        atom_index_sets = ((picked, signs)
                           for size in range(1, atom_degree+1)
                           for picked in combinations(range(dimension), size)
                           for signs in product((0, 1), repeat=size))
    else:
        literals = [(coordinate, sign)
                    for coordinate in range(dimension) for sign in (0, 1)]
        atom_index_sets = ((tuple(literals[t][0] for t in picked), tuple(literals[t][1] for t in picked))
                           for size in (1, 2)
                           for picked in combinations_with_replacement(range(len(literals)), size))
    for picked, signs in atom_index_sets:
        atom = const(1, dimension)
        for coordinate, sign in zip(picked, signs):
            factor = variable(coordinate, dimension)
            if sign:
                factor = add(const(1, dimension), factor, Fraction(-1))
            atom = mul(atom, factor)
        value_of_atom = evaluate(functional, atom)
        if value_of_atom < 1:
            raise SystemExit('functional is below one on an allowed positive atom')
        atom_values.append(value_of_atom)
    print(json.dumps({
        'status': 'verified',
        'annihilated_multiplier_columns': annihilated,
        'positive_atoms_checked': len(atom_values),
        'minimum_atom_value': [min(atom_values).numerator, min(atom_values).denominator],
        'conclusion': 'No nonzero nonnegative atom combination equals a linear combination of the declared multiplier columns.'
    }, indent=2))


if __name__ == '__main__':
    main()
