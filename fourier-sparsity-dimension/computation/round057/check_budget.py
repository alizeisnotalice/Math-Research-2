"""Exact checks of rank, slice occupancy and the field example (standard library)."""
from fractions import Fraction
from collections import Counter
from itertools import product
import json
import platform


def walsh(values):
    values = list(values)
    stride = 1
    while stride < len(values):
        for start in range(0, len(values), 2 * stride):
            for i in range(start, start + stride):
                a, b = values[i], values[i + stride]
                values[i], values[i + stride] = a+b, a-b
        stride *= 2
    return values


def support(values):
    return [i for i, x in enumerate(walsh(values)) if x]


def rank(points):
    pivots = {}
    for value in points:
        while value:
            j = value.bit_length()
            if j not in pivots:
                pivots[j] = value
                break
            value ^= pivots[j]
    return len(pivots)


def audit(table, inner):
    rows = [support(table[i:i+inner]) for i in range(0,len(table),inner)]
    counts = Counter(b for row in rows for b in row)
    assert all(len(row) == 2**(rank([inner+b for b in row])-1) for row in rows)
    sparsity = len(support(table))
    occ = sum(Fraction(len(rows),v) for v in counts.values())
    projected = rank([inner+b for b in counts])
    assert 4*projected**2 <= 9*occ <= 9*sparsity
    return dict(s=sparsity, lifted_rank=projected, occupancy=str(occ))


def field_prod(x,y):
    z = 0
    while y:
        if y & 1:
            z ^= x
        y >>= 1
        x <<= 1
        if x & 4:
            x ^= 7
    return z


def main():
    tight = 0
    for mask in range(65536):
        v=[1-2*((mask >> j)&1) for j in range(16)]
        a=audit(v,4)
        tight += 4*a['lifted_rank']**2 == 9*a['s']
    field = [1-2*(field_prod(x%4,y%4)==field_prod(x//4,y//4))
             for y in range(16) for x in range(16)]
    example = audit(field,16)
    assert example['s']==256 and rank(support(field))==8
    models = []
    for m in (1,2,4):
        inner=2**(3*m)
        v = [1-2*(((x>>(3*y))&1)^(((x>>(3*y+1))&3)==3))
             for y in range(m) for x in range(inner)]
        a=audit(v,inner)
        assert a['lifted_rank']==3*m and a['s']==4*m*m
        models.append(a)
    print(json.dumps(dict(python=platform.python_version(), four_variable_tables=65536,
                          equality_count=tight, field_q4_m2=example,
                          private_plane_models=models),indent=2))


if __name__ == '__main__':
    main()
