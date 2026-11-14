"""Exact finite receiver LP checks; Python standard library only.

Run: python3 computation/general_dual_obstructions.py
No numerical tolerances or random input; all arithmetic is Fraction.
"""
from fractions import Fraction as F
from itertools import combinations, product
import json
import platform


def solve(a, b):
    n = len(b)
    a = [list(map(F, row)) + [F(v)] for row, v in zip(a, b)]
    for k in range(n):
        p = next((j for j in range(k, n) if a[j][k]), None)
        if p is None:
            return None
        a[k], a[p] = a[p], a[k]
        q = a[k][k]
        a[k] = [v / q for v in a[k]]
        for j in range(n):
            if j != k:
                q = a[j][k]
                a[j] = [v - q*w for v, w in zip(a[j], a[k])]
    return tuple(row[-1] for row in a)


def lp(masks, objective):
    n = len(objective)
    rows = [tuple(F(i in mask) for i in range(n)) for mask in masks]
    rhs = [F(1)] * len(rows)
    rows += [tuple(F(-int(i == j)) for i in range(n)) for j in range(n)]
    rhs += [F(0)] * n
    best, arg = F(-1), None
    for active in combinations(range(len(rows)), n):
        z = solve([rows[i] for i in active], [rhs[i] for i in active])
        if z is None or any(sum(a*b for a,b in zip(row,z)) > r for row,r in zip(rows,rhs)):
            continue
        value = sum(a*b for a,b in zip(objective,z))
        if value > best:
            best, arg = value, z
    return best, arg


def contains(q, p):
    x,y,r = q
    return abs(x-p[0]) <= r and abs(y-p[1]) <= r


def global_masks(squares):
    axes = []
    for dim in (0,1):
        ends = sorted({q[dim]+sgn*q[2] for q in squares for sgn in (-1,1)})
        axes.append(ends + [(a+b)/2 for a,b in zip(ends,ends[1:])])
    masks = {frozenset(i for i,q in enumerate(squares) if contains(q,p)) for p in product(*axes)}
    # Nonnegative receiver weights make every proper-subset constraint redundant.
    return sorted((s for s in masks if s and not any(s < t for t in masks)), key=lambda s: tuple(sorted(s)))


def main():
    n = 5
    squares = [(F(i,n+1),1-F(i,n+1),F(1)) for i in range(1,n+1)]
    sources = [(x+1,y+1) for x,y,r in squares]
    source_masks = [frozenset(i for i,q in enumerate(squares) if contains(q,p)) for p in sources]
    assert source_masks == [frozenset([i]) for i in range(n)]
    assert all(contains(q,(F(0),F(0))) for q in squares)
    support_value, support_z = lp(source_masks,[F(4)]*n)
    all_masks = global_masks(squares)
    all_value, all_z = lp(all_masks,[F(4)]*n)
    assert support_value == 4*n and all_value == 4
    assert all(sum(F(5) for p in sources if contains(q,p)) > 4 for q in squares)

    nested = [(F(0),F(0),F(1))] + [(F(0),F(0),F(1,n)) for _ in range(n-1)]
    areas = [4*q[2]**2 for q in nested]
    nested_value, nested_z = lp(global_masks(nested), areas)
    equal_weight_value = sum(areas)/n
    assert nested_value == 4 and equal_weight_value == F(116,125)
    # Two unit atoms at (+/-1/2,0): full halo = singleton union + pair protrusions.
    # Exact quadratic-field pairs represent a+b*sqrt(2).
    def mul(a,b):
        return (a[0]*b[0]+2*a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    extra = mul((-1,1),(-1,1))
    assert extra == (3,-2) and extra[0]**2 > 2*extra[1]**2
    result = {
        'python':platform.python_version(), 'arithmetic':'fractions.Fraction; exhaustive LP vertices',
        'seed':None, 'tolerance':None,
        'support_gap': {'n':n, 'squares':squares, 'sources':sources, 'source_mass_each':5,
                        'support_value':support_value, 'support_z':support_z,
                        'all_space_value':all_value,'all_space_z':all_z,
                        'maximal_global_masks':[sorted(s) for s in all_masks]},
        'unweighted_nerve_rule':{'n':n,'areas':areas,'exact_optimum':nested_value,
                                'equal_weight_value':equal_weight_value,
                                'loss_factor':nested_value/equal_weight_value},
        'singleton_halo_loss':{'source_points':[[-F(1,2),0],[F(1,2),0]],
                               'masses':[1,1], 'singleton_union_area':2,
                               'extra_area_a_plus_b_sqrt2':extra,
                               'complete_halo_area_a_plus_b_sqrt2':[5,-2]},
        'all_assertions_passed':True,
    }
    print(json.dumps(result,default=str,indent=2))


if __name__ == '__main__':
    main()
