#!/usr/bin/env python3
"""Exact fixtures where all three exterior components occur at once."""
from collections import Counter
import json
import platform
import time

from check_components import audit_support, compute_components, gaussian, points, walsh_mask


def sign(value):
    return 1 - 2 * value


def fixture(bent_count):
    private = 2*bent_count
    n = private+2
    core = (1 << (1 << private))-1
    values = []
    for x in range(1<<n):
        phases = [sign(((x >> (2*j)) & 1)*((x >> (2*j+1)) & 1)) if j<bent_count else 1
                  for j in range(3)]
        a,b,c = phases
        y1,y2 = ((x>>private)&1),((x>>(private+1))&1)
        numerator = -a*b*c+a*sign(y1)+b*sign(y2)+c*sign(y1^y2)
        if numerator not in (-2,2):
            raise AssertionError('fixture not Boolean')
        values.append(numerator//2)
    truth = sum((w==-1)<<x for x,w in enumerate(values))
    support = walsh_mask(n,truth)
    direct = [sum(w*sign(bin(x&u).count('1')%2) for x,w in enumerate(values)) for u in range(1<<n)]
    support_direct = sum((a!=0)<<u for u,a in enumerate(direct))
    if support != support_direct:
        raise AssertionError('direct Walsh and butterfly supports differ')
    expected_total = 4**bent_count+4*bent_count+(3-bent_count)
    if len(tuple(points(support))) != expected_total:
        raise AssertionError('unexpected full support')
    blocks = sorted((tuple(points(c)) for c in compute_components(n,support&~core)),key=lambda row: row[0]>>private)
    projection = [sorted({x>>private for x in row}) for row in blocks]
    expected_sizes = [4 if j<bent_count else 1 for j in range(3)]
    if len(blocks)!=3 or projection != [[1],[2],[3]] or [len(row) for row in blocks] != expected_sizes:
        raise AssertionError('incorrect triple block classification')
    census = None
    if n<=6:
        counters,_,_ = audit_support(n,support)
        census = [dict(codimension=t,components=k,count=count) for (t,k),count in sorted(counters.items())]
        if sum(r['count'] for r in census) != sum(gaussian(n,i) for i in range(n+1)):
            raise AssertionError('incomplete fixture subspace audit')
    return dict(n=n,dimension_of_core=private,truth_mask=truth,sparsity=expected_total,
                component_frequencies=blocks,quotient_labels=[1,2,3],
                full_walsh_nonzero_coefficients=[[u,direct[u]] for u in points(support)],
                all_subspace_census=census)


def main():
    start = time.perf_counter()
    records = [fixture(j) for j in range(4)]
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),
                          records=records,wall_seconds=time.perf_counter()-start),indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
