"""Exact Boolean model for the four active common labels in the proof."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
import json
import platform


def walsh(row):
    row = row.copy()
    k = 1
    while k < len(row):
        for i in range(0, len(row), 2*k):
            for j in range(i, i+k):
                a, b = row[j], row[j+k]
                row[j], row[j+k] = a+b, a-b
        k *= 2
    return row


def dim(vectors):
    pivots = {}
    for x in vectors:
        while x and x.bit_length() in pivots:
            x ^= pivots[x.bit_length()]
        if x:
            pivots[x.bit_length()] = x
    return len(pivots)


def model(v, z):
    target = 4 - 2*v.bit_count()
    if abs(target) == 4:
        return 1 if target > 0 else -1
    if target:
        return (1 if target > 0 else -1) * (-1 if z & 3 == 3 else 1)
    return (-1 if z & 12 == 12 else 1) * (-1 if z & 16 else 1)


if __name__ == '__main__':
    rows = [[model(v,z) for z in range(32)] for v in range(16)]
    mean = [Q(sum(row), 32) for row in rows]
    costs = {Q(1):Q(0), Q(1,2):Q(4,3), Q(0):Q(9,4)}
    joint = 4 + sum(costs[abs(g)] for g in mean)/16
    assert joint == Q(529,96)
    common = [b for b, x in enumerate(walsh(mean)) if x]
    assert common == [1,2,4,8]
    labels = {b for row in rows for b, x in enumerate(walsh(row)) if x}
    mult = Counter(a^b for a,b in combinations(labels,2))
    planes = [v for v in combinations(sorted(labels),4) if v[0]^v[1]^v[2]^v[3] == 0]
    assert all(min(mult[v[0]^b] for b in v[1:]) == 2 for v in planes)
    spectrum = [i for i,x in enumerate(walsh([x for row in rows for x in row])) if x]
    rank, s = dim(spectrum), len(spectrum)
    assert 96*rank**2 <= 529*s
    print(json.dumps(dict(python=platform.python_version(), variables=9, common=common,
                          projected_labels=sorted(labels), planes=planes,
                          sparsity=s, rank=rank, joint_local_cost=str(joint)),indent=2))
