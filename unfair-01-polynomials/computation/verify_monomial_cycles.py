#!/usr/bin/env python3
"""Independent rational verifier. Does not import the searcher or use floats."""
from fractions import Fraction
import argparse
import json


def verify(record):
    I, J = record['I'], record['J']
    for S in (I, J):
        assert S and S[0] == 0 and S == sorted(set(S))
        assert all(type(i) is int and i >= 0 for i in S)
    def edges(k):
        return [(i, k-i) for i in I if k-i in J]
    unit = {('a', 0), ('a', I[-1]), ('b', 0), ('b', J[-1])}
    for k in range(I[-1]+J[-1]+1):
        es = edges(k)
        if len(es) == 1:
            i, j = es[0]
            unit.update([('a', i), ('b', j)])
    variables = sorted([('a', i) for i in I if ('a', i) not in unit] +
                       [('b', j) for j in J if ('b', j) not in unit])
    # Every remaining variable is <1 by its own positive-index budget.
    assert all(len(edges(v[1])) > 1 for v in variables)
    def monomials(k):
        return sorted(tuple(v for v in [('a', i), ('b', j)] if v not in unit)
                      for i, j in edges(k))
    total = {v: Fraction() for v in variables}
    assert record['terms']
    for term in record['terms']:
        k, ell = term['diagonals']
        assert type(k) is int and type(ell) is int and k != ell
        left, right = monomials(k), monomials(ell)
        assert left and right
        # List removal preserves multiplicity independently of Counter subtraction.
        for monomial in left[:]:
            if monomial in right:
                left.remove(monomial)
                right.remove(monomial)
        assert len(left) == len(right) == 1
        assert isinstance(term['weight'], str), 'Use exact rational strings'
        weight = Fraction(term['weight'])
        for v in left[0]:
            total[v] += weight
        for v in right[0]:
            total[v] -= weight
    assert all(t >= 0 for t in total.values()) and any(t > 0 for t in total.values())
    assert [list(v) for v in variables] == record['variables']
    assert [str(total[v]) for v in variables] == record['z']
    return True


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('certificates')
    args = ap.parse_args()
    with open(args.certificates) as handle:
        records = json.load(handle)
    checked = 0
    for record in records:
        if 'terms' in record:
            verify(record)
            checked += 1
    print(json.dumps({'independently_verified': checked, 'records': len(records)}))
