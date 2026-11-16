#!/usr/bin/env python3
"""Exact prescribed-peak-pair audit for D4/reflection candidates.

The C++ kernel is the project pair-DP, compiled into this temporary directory.
It counts words with exactly the specified two (nonadjacent) peaks.  This
script uses only its integer stdout and performs no floating-point decisions.
"""

from collections import Counter
from fractions import Fraction
import itertools
import json
import subprocess
import sys


def xy(n, v):
    return (v % n, v // n)


def d1(n, p, q):
    x, y = xy(n, p); u, v = xy(n, q)
    return abs(x-u) + abs(y-v)


def refl_x(n, v):
    x, y = xy(n, v)
    return y*n + (n-1-x)


def pair_count(binary, n, p, q):
    return int(subprocess.check_output(
        [binary, str(n), str(p), str(q)], text=True))


def audit(binary, n):
    N = n*n
    K = {}
    for p, q in itertools.combinations(range(N), 2):
        if d1(n, p, q) > 1:
            K[p, q] = pair_count(binary, n, p, q)

    def k(p, q):
        return K[tuple(sorted((p, q)))]

    total = sum(K.values())
    expected = {3: 143808, 4: 779299174080}[n]
    assert total == expected
    cx = Fraction(n-1, 2)
    vx = vy = cxx = cyy = cxy = Fraction(0)
    # K is unordered; use the symmetric ordered-pair law in the products.
    for (p, q), value in K.items():
        xp, yp = xy(n, p); xq, yq = xy(n, q)
        vx += value*((xp-cx)**2 + (xq-cx)**2)
        vy += value*((yp-cx)**2 + (yq-cx)**2)
        cxx += 2*value*(xp-cx)*(xq-cx)
        cyy += 2*value*(yp-cx)*(yq-cx)
        cxy += value*((xp-cx)*(yq-cx) + (xq-cx)*(yp-cx))
    vx /= 2*total; vy /= 2*total
    cxx /= 2*total; cyy /= 2*total; cxy /= 2*total

    # Orbit-summed version of the four-point reflection comparison.  It is
    # deliberately an aggregate statement: pointwise comparison is tested
    # separately below and fails.
    same_side = cross_side = 0
    for (p, q), value in K.items():
        xp, _ = xy(n, p); xq, _ = xy(n, q)
        sp = (xp > cx) - (xp < cx)
        sq = (xq > cx) - (xq < cx)
        if sp and sq:
            if sp == sq:
                same_side += value
            else:
                cross_side += value

    # Stronger than the preceding total comparison: retain the left-hand root
    # p and sum only over q in the same strict half-plane, comparing q with
    # its reflected mate.  This is still an aggregate four-point statement,
    # but need not follow from D4.
    per_anchor = Counter()
    anchor_examples = {}
    for p in range(N):
        xp, _ = xy(n, p)
        if xp >= cx:
            continue
        same = cross = 0
        for q in range(N):
            xq, _ = xy(n, q)
            if xq >= cx or q == p or d1(n, p, q) <= 1:
                continue
            r = refl_x(n, q)
            if r == p or d1(n, p, r) <= 1:
                continue
            same += k(p, q)
            cross += k(p, r)
        if same or cross:
            verdict = 'holds' if cross >= same else 'fails'
            per_anchor[verdict] += 1
            anchor_examples.setdefault(verdict, {
                'anchor_p': list(xy(n, p)), 'same_side_sum': same,
                'reflected_cross_sum': cross})

    # Candidate: reflecting q across the central vertical axis toward p should
    # not lower the exact two-peak kernel.  Check every eligible ordered pair.
    comparisons = Counter()
    examples = {}
    for p, q in K:
        for a, b in ((p, q), (q, p)):
            r = refl_x(n, b)
            if r == a or d1(n, a, r) <= 1:
                continue
            if d1(n, a, r) >= d1(n, a, b):
                continue
            left, right = k(a, r), k(a, b)
            verdict = 'holds' if left >= right else 'fails'
            comparisons[verdict] += 1
            examples.setdefault(verdict, {
                'p': list(xy(n, a)), 'far_q': list(xy(n, b)),
                'reflected_q': list(xy(n, r)),
                'far_distance': d1(n, a, b),
                'near_distance': d1(n, a, r),
                'K_near': left, 'K_far': right,
            })
    assert total > 0
    return {
        'n': n, 'A2': total,
        'ordered_root_covariance_about_center': {
            'Var_Xx': str(vx), 'Var_Xy': str(vy),
            'Cxx': str(cxx), 'Cyy': str(cyy), 'Cxy': str(cxy),
            'dot_covariance': str(cxx + cyy),
            'E_norm_X_minus_c_sq': str(vx + vy),
            'E_root_distance_sq_from_variance_identity':
                str(2*(vx + vy - cxx - cyy)),
            'Q_times_n_sq': str(Fraction(sum(
                value*((xy(n, p)[0]-xy(n, q)[0])**2 +
                       (xy(n, p)[1]-xy(n, q)[1])**2)
                for (p, q), value in K.items()), total)),
        },
        'orbit_summed_vertical_reflection': {
            'definition': 'strict left/right only; compare total cross-side peak-pair mass to same-side mass',
            'same_side_mass': same_side, 'cross_side_mass': cross_side,
            'cross_minus_same': cross_side - same_side,
            'holds_cross_ge_same': cross_side >= same_side,
        },
        'anchor_summed_vertical_reflection': {
            'definition': 'for each strict-left anchor p, sum eligible q in strict-left and compare K(p,refl(q)) to K(p,q)',
            'anchors': dict(per_anchor), 'examples': anchor_examples,
        },
        'toward_reflection_candidate': {
            'definition': 'if x-reflection of q strictly reduces d1(p,q), require K(p,refl(q)) >= K(p,q)',
            'comparisons': dict(comparisons), 'examples': examples},
    }


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: audit_pair_reflection.py PAIR_DP_BINARY OUTPUT_JSON')
    result = {'source': 'exact prescribed peak-pair DP',
              'results': [audit(sys.argv[1], n) for n in (3, 4)]}
    with open(sys.argv[2], 'w') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps(result, indent=2))
