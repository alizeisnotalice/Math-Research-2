#!/usr/bin/env python3
"""Propose exact row-space product certificates; LP failure stays unknown."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
import argparse
import json
import numpy as np
from scipy.optimize import linprog


def rows_for_support(I, J):
    edges = defaultdict(list)
    for i in I:
        for j in J:
            edges[i+j].append((('a', i), ('b', j)))
    units = {('a', 0), ('a', I[-1]), ('b', 0), ('b', J[-1])}
    for es in edges.values():
        if len(es) == 1:
            units.update(es[0])
    variables = sorted(({('a', i) for i in I} | {('b', j) for j in J}) - units)
    reduced = {k: Counter(tuple(int(v in e) for v in variables) for e in es)
               for k, es in edges.items()}
    rows = []
    for k, ell in combinations(sorted(edges), 2):
        left, right = reduced[k]-reduced[ell], reduced[ell]-reduced[k]
        if sum(left.values()) == sum(right.values()) == 1:
            a, b = next(iter(left)), next(iter(right))
            rows.append({'diagonals': [k, ell], 'row': [x-y for x, y in zip(a, b)]})
    return variables, rows


def find_certificate(I, J):
    variables, rows = rows_for_support(I, J)
    result = {'I': I, 'J': J, 'status': 'unknown'}
    if not rows or not variables:
        return result
    L = np.array([r['row'] for r in rows], dtype=float)
    # Normalize the nonzero nonnegative vector yL by its coordinate sum.
    fit = linprog(np.zeros(len(rows)), A_ub=-L.T, b_ub=np.zeros(len(variables)),
                  A_eq=L.sum(axis=1)[None, :], b_eq=[1],
                  bounds=[(None, None)]*len(rows), method='highs')
    result['lp_status'] = int(fit.status)
    if not fit.success:
        return result
    for bound in (100, 10000, 1000000):
        y = [Fraction(float(v)).limit_denominator(bound) for v in fit.x]
        z = [sum((y[i]*rows[i]['row'][j] for i in range(len(rows))), Fraction())
             for j in range(len(variables))]
        if any(z) and all(v >= 0 for v in z):
            result.update(status='proposed_product_certificate',
                          terms=[{'diagonals': row['diagonals'], 'weight': str(w)}
                                 for row, w in zip(rows, y) if w],
                          variables=[list(v) for v in variables], z=[str(v) for v in z])
            return result
    result['reason'] = 'floating proposal did not reconstruct an exact certificate'
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('input', help='JSON array of {I,J} support records')
    args = ap.parse_args()
    with open(args.input) as handle:
        cases = json.load(handle)
    print(json.dumps([find_certificate(c['I'], c['J']) for c in cases], indent=2))
