"""Exact checks of parity-circuit laws and two finite-valued Boolean lifts."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations
from math import comb
import json
import platform


def walsh(values):
    a = list(values)
    width = 1
    while width < len(a):
        for start in range(0, len(a), 2*width):
            for i in range(start, start+width):
                x, y = a[i], a[i+width]
                a[i], a[i+width] = x+y, x-y
        width *= 2
    return a


def laws():
    rows = []
    for budgets in [(1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,)]:
        distribution = {Q(0): Q(1)}
        common = Q(0)
        for b in budgets:
            m = 2*b+2
            new = defaultdict(Q)
            for old, p in distribution.items():
                for j in range(1, m, 2):
                    new[old+Q(m-2*j, 4)] += p*Q(comb(m, j), 2**(m-1))
            distribution = new
            common += Q((m-1)**2, m)
        private = Q(9, 4)*distribution[0] + Q(4, 3)*(distribution[1]+distribution[-1])
        rows.append({"budgets": budgets, "common": str(common), "private": str(private),
                     "joint": str(common+private)})
    assert max(Q(x["joint"]) for x in rows) == Q(1009, 96)
    return rows


def actual_lift(common_dim, private_dim, common_labels, row_builder):
    rows = [row_builder(x) for x in range(1 << common_dim)]
    zero_band = [sum(row) for row in rows]
    w = walsh(zero_band)
    assert {j for j, a in enumerate(w) if a} == set(common_labels)
    qset = set()
    norms = Counter()
    for row in rows:
        assert len(row) == 1 << private_dim and set(row) <= {-1, 1}
        spectrum = walsh(row)
        qset.update(j for j, a in enumerate(spectrum) if a)
        norms[str(abs(Q(sum(row), len(row))))] += 1
    pairs = Counter(u ^ v for u, v in combinations(qset, 2))
    assert max(pairs.values()) <= 2
    return {"physical_rows": len(rows), "truth_table_size": len(rows)*len(rows[0]),
            "common_size": len(common_labels), "abs_mean_histogram": norms,
            "private_labels": sorted(qset),
            "max_private_unordered_difference_multiplicity": max(pairs.values())}


def sharp_ten():
    tags = [1 << i for i in range(9)] + [511]
    def row(x):
        numerator = sum((-1)**((x & v).bit_count() % 2) for v in tags[:-1])
        numerator -= (-1)**(x.bit_count() % 2)
        g = Q(numerator, 8)
        negatives = 2*(1-g)
        assert negatives.denominator == 1
        return [-1 if i < negatives else 1 for i in range(4)]
    return actual_lift(9, 2, tags, row)


def sharp_sixteen():
    groups = [[j << (3*i) for j in (1, 2, 4, 7)] for i in range(4)]
    tags = sum(groups, [])
    def row(x):
        total = 0
        for group in groups:
            signs = [(-1)**((x & u).bit_count() % 2) for u in group]
            total += sum(signs[:3]) - signs[3]
        g = Q(total, 8)
        result = []
        for z in range(32):
            if abs(g) == 1:
                val = int(g)
            elif abs(g) == Q(1, 2):
                val = int(2*g) * (-1)**int((z & 3) == 3)
            else:
                assert g == 0
                val = (-1)**(((z >> 2) & 1) ^ int((z & 24) == 24))
            result.append(val)
        return result
    return actual_lift(12, 5, tags, row)


print(json.dumps({"python": platform.python_version(), "exact_joint_costs": laws(),
                  "ten_label_conditional_mean": sharp_ten(),
                  "sixteen_label_conditional_mean": sharp_sixteen()}, indent=2))
