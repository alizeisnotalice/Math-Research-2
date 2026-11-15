#!/usr/bin/env python3
"""Independent Fraction polynomial arithmetic and Sturm certificate verifier."""
from fractions import Fraction as Q
import argparse
import json


def trim(p):
    p = list(map(Q, p))
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(a, b):
    return trim([(a[i] if i < len(a) else 0)+(b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def scale(a, c):
    return trim([x*c for x in a])


def mul(a, b):
    p = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            p[i+j] += x*y
    return trim(p)


def remainder(a, b):
    p = trim(a)
    while p != [0] and len(p) >= len(b):
        shift, c = len(p)-len(b), p[-1]/b[-1]
        for j, x in enumerate(b):
            p[j+shift] -= c*x
        p = trim(p)
    return p


def evaluate(p, x):
    y = Q(0)
    for a in reversed(p):
        y = y*x+a
    return y


def sturm(p):
    p = trim(p)
    seq = [p, trim([i*p[i] for i in range(1, len(p))])]
    while seq[-1] != [0]:
        nxt = scale(remainder(seq[-2], seq[-1]), -1)
        if nxt == [0]:
            break
        seq.append(scale(nxt, 1/abs(nxt[-1])))
    return seq


def variations(seq, x):
    signs = [(v > 0)-(v < 0) for p in seq if (v := evaluate(p, x)) != 0]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def verify(case):
    r = case['r']
    assert type(r) is int and r in (3, 5, 7, 9, 11)
    t = [Q((-1)**i) for i in range(r+1)]
    b = [[Q(1)]]
    for k in range(1, 2*r):
        p = add([1], scale(mul([0, 1], b[-1]), -1))
        if k >= r:
            p = add(p, scale(mul(t, b[k-r]), -1))
        b.append(p)
    F = b[-1]
    H = trim(case['H'])
    assert mul([0, -1, 1], H) == F
    R = add(mul([-1, 2], add([0, 1], scale(t, -1))),
            scale(mul(mul(t, t), [1, -3, 3]), -1))
    assert R == trim(case['R'])
    lo, hi = map(Q, case['interval'])
    assert 0 < lo < hi < 1
    assert evaluate(H, 0)*evaluate(H, 1) != 0
    assert evaluate(H, lo)*evaluate(H, hi)*evaluate(R, lo)*evaluate(R, hi) != 0
    hs, rs = sturm(H), sturm(R)
    hv = [variations(hs, x) for x in [Q(0), Q(1), lo, hi]]
    rv = [variations(rs, x) for x in [lo, hi]]
    assert hv[0]-hv[1] == hv[2]-hv[3] == 1
    assert rv[0]-rv[1] == 0
    assert evaluate(R, (lo+hi)/2) > 0
    return {'r': r, 'H_variations_at_0_1_lo_hi': hv,
            'R_variations_at_lo_hi': rv, 'R_midpoint_sign': 1}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('certificate')
    args = ap.parse_args()
    with open(args.certificate) as handle:
        data = json.load(handle)
    assert [c['r'] for c in data['small_odd']] == [3, 5, 7, 9, 11]
    print(json.dumps([verify(c) for c in data['small_odd']], indent=2))
