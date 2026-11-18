#!/usr/bin/env python3
"""Exact degree-eight S={3} coefficient bounds and prefix replay."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import isqrt
from pathlib import Path

ZERO, ONE = (0, 0), (2, 0)
BOUND = {3: (65, 2), 4: (257, 4), 5: (1025, 8),
         6: (4097, 16), 7: (16385, 32), 8: (65537, 64)}
BOX = {4: 120, 5: 92, 6: 41, 7: 9, 8: 1}


def ci(x):
    return 2*x, 0


def minus(x):
    return -x[0], -x[1]


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def mult(x, y, D):
    u, v = x
    p, q = y
    r, s = u*p+v*q*D, u*q+v*p
    assert r % 2 == s % 2 == 0
    return r//2, s//2


def scale(n, x):
    return n*x[0], n*x[1]


def conjugate(x):
    return x[0], -x[1]


def sign(x, D):
    u, v = x
    if not v:
        return (u > 0)-(u < 0)
    if not u or u*v > 0:
        return (v > 0)-(v < 0)
    difference = u*u-v*v*D
    assert difference != 0
    return ((difference > 0)-(difference < 0))*((u > 0)-(u < 0))


def shifts(x, D):
    u, v = x
    if not v:
        assert u % 2 == 0
        return ((-u//2, ZERO), (1-u//2, ONE))
    radical_floor = isqrt(v*v*D)
    if v < 0:
        radical_floor = -radical_floor-1
    shift = -((u+radical_floor)//2)
    shifted = add(x, ci(shift))
    assert sign(shifted, D) > 0 and sign(add(ONE, minus(shifted)), D) > 0
    return ((shift, shifted),)


def choose(base, row, D):
    if row == 2:
        return ((0, base),) if base in (ZERO, ONE) else ()
    if row == 5:
        return ((-base[0]//2, ZERO),) if not base[1] and base[0] % 2 == 0 else ()
    return shifts(base, D)


def coefficient_without_current(poly, prefix, degree):
    value = ZERO
    for j in range(1, degree+1):
        value = add(value, scale(prefix[degree-j], poly[j]))
    return value


INITIAL = (((1,), (1,), (ONE,), (ONE,)),)


def row_step(states, poly, degree, D):
    pb = tuple(map(conjugate, poly))
    valid = []
    for G, Q, A, B in states:
        a = coefficient_without_current(poly, G, degree)
        b = coefficient_without_current(pb, Q, degree)
        cross = ZERO
        for i in range(1, degree):
            cross = add(cross, mult(A[i], B[degree-i], D))
        for g, ak in choose(a, degree, D):
            for q, bk in choose(b, degree, D):
                if add(cross, add(ak, bk)) in (ZERO, ONE):
                    valid.append((G+(g,), Q+(q,), A+(ak,), B+(bk,)))
    return valid


def strict_house(x, num, den, D):
    u, v = x
    gap = 2*num-den*abs(u)
    return gap > 0 and gap*gap > den*den*v*v*D


def newton_without_current(poly, powers, degree, D):
    term = ZERO
    for j in range(1, degree):
        term = add(term, mult(poly[j], powers[degree-j], D))
    return minus(term)


def newton_next(poly, powers, degree, D):
    return add(newton_without_current(poly, powers, degree, D), scale(-degree, poly[degree]))


def integer_interval(base, bound, slope, D):
    num, den = bound
    u, v = base
    radius = (2*num-isqrt((den*v)**2*D)-1)//den
    if radius < 0:
        return 1, 0
    assert slope < 0
    step = -2*slope
    return -((radius-u)//step), (u+radius)//step


def strict_integer_ball(base, degree, D):
    lo, hi = integer_interval(base, BOUND[degree], -degree, D)
    return range(max(-BOX[degree], lo), min(BOX[degree], hi)+1)


def prepare_translated_row(states, poly, degree, D):
    """Evaluate with the integer coefficient of x^(degree-3) set to zero."""
    pb = tuple(map(conjugate, poly))
    prepared = []
    L, M = poly[1][0]//2, poly[2][0]//2
    assert M == L*L
    for G, Q, A, B in states:
        cross = ZERO
        for i in range(1, degree):
            cross = add(cross, mult(A[i], B[degree-i], D))
        left = coefficient_without_current(poly, G, degree)
        right = coefficient_without_current(pb, Q, degree)
        slope_a = add(ci(G[3]+L*M), minus(poly[3]))
        slope_b = add(ci(Q[3]+L*M), minus(pb[3]))
        prepared.append((cross, left, right, slope_a, slope_b))
    return prepared


def count_translated_row(prepared, coefficient, D):
    count = 0
    for cross, left, right, slope_a, slope_b in prepared:
        for _, ak in shifts(add(left, scale(coefficient, slope_a)), D):
            for _, bk in shifts(add(right, scale(coefficient, slope_b)), D):
                count += add(cross, add(ak, bk)) in (ZERO, ONE)
    return count


def ramp_sum(c, step, lo, hi):
    """Sum max(0, c+step*i) over an integer interval."""
    lo = max(lo, (-c)//step+1)
    if lo > hi:
        return 0
    n = hi-lo+1
    return c*n + step*n*(lo+hi)//2


def tail_box(poly, powers, D):
    """Count admissible (d,e), grouping e intervals by the slope of s7."""
    L = poly[1][0]//2
    assert poly[2] == ci(L*L) and poly[6] == poly[7] == ZERO
    powers = list(powers)
    powers.append(newton_next(poly, powers, 6, D))
    powers.append(newton_next(poly, powers, 7, D))
    p8 = newton_next(poly, powers, 8, D)
    dlo, dhi = integer_interval(powers[6], BOUND[6], -6, D)
    dlo, dhi = max(-41,dlo), min(41,dhi)
    zlo, zhi = integer_interval(powers[7], BOUND[7], -7, D)
    if dlo > dhi or zlo > zhi:
        return None
    if L == 0:
        if not strict_house(p8, *BOUND[8], D):
            return None
        elo, ehi = max(-9,zlo), min(9,zhi)
        return ((dhi-dlo+1)*(ehi-elo+1), dlo, dhi, elo, ehi) if elo <= ehi else None
    slope = abs(L)
    elo, ehi = integer_interval(p8, BOUND[8], -8*slope, D)
    if L > 0:
        elo, ehi = -ehi,-elo
    elo, ehi = max(-9,elo), min(9,ehi)
    if elo > ehi:
        return None
    if L < 0:
        dlo, dhi = -dhi,-dlo
    dlo = max(dlo, (elo-zhi+slope-1)//slope)
    dhi = min(dhi, (ehi-zlo)//slope)
    if dlo > dhi:
        return None
    count = (ramp_sum(zhi-elo+1,slope,dlo,dhi) - ramp_sum(zhi-ehi,slope,dlo,dhi)
             - ramp_sum(zlo-elo,slope,dlo,dhi) + ramp_sum(zlo-ehi-1,slope,dlo,dhi))
    assert count > 0
    support_d = (dlo,dhi) if L > 0 else (-dhi,-dlo)
    return count, *support_d, max(elo,zlo+slope*dlo), min(ehi,zhi+slope*dhi)


def cores():
    pairs = [(L,L*L) for L in range(-4,5)]
    entries = []
    for L, M in pairs:
        for T in range(-184,185):
            gap = 65-abs(-2*L**3+6*L*M-3*T)
            if gap <= 0:
                continue
            largest_D = min((gap*gap-1)//9, (185-abs(T))**2-1)
            for D in range(T*T % 4, largest_D+1, 4):
                if D <= 0 or isqrt(D)**2 == D:
                    continue
                N = (T*T-D)//4
                entries.append((L,M,T,N,D))
    return entries


def audit_core(job):
    core, max_row, final = job
    L, M, T, N, D = core
    poly = [ONE, ci(L), ci(M), (T, 1), ZERO, ZERO, ZERO, ZERO, ONE]
    powers = [ci(8)]
    for k in range(1, 4):
        powers.append(newton_next(poly, powers, k, D))
    assert strict_house(poly[3], 185, 2, D)
    assert strict_house(powers[2], 17, 1, D)
    assert strict_house(powers[3], 65, 2, D)
    states = INITIAL
    for k in range(1,4):
        states = row_step(states, poly, k, D)
    assert len(states) == (3 if L in (0,1) else 1)
    counts = Counter(rows3=len(states))
    ranges = {"core_L": [L,L], "core_M": [M,M], "core_T": [T,T],
              "core_N": [N,N], "core_D": [D,D]}
    sha = hashlib.sha256()

    def record(name, value):
        interval = ranges.setdefault(name, [value,value])
        interval[0] = min(interval[0], value)
        interval[1] = max(interval[1], value)

    base4 = newton_without_current(poly, powers, 4, D)
    for b in strict_integer_ball(base4, 4, D):
        poly[4] = ci(b)
        p4 = add(base4, ci(-4*b))
        counts["newton4_b"] += 1
        s4 = row_step(states, poly, 4, D)
        counts["rows4"] += len(s4)
        if not s4:
            continue
        s5 = row_step(s4, poly, 5, D)
        s6 = row_step(s5, poly, 6, D)
        s7 = row_step(s6, poly, 7, D)
        counts.update(rows5=len(s5), rows6=len(s6), rows7=len(s7))
        if not s7:
            continue
        counts["row7_core_b"] += 1
        record("row7_b", b)
        if max_row == 7:
            continue
        base5 = newton_without_current(poly, [*powers,p4], 5, D)
        prepared = prepare_translated_row(s7, poly, 8, D)
        for c in strict_integer_ball(base5, 5, D):
            poly[5] = ci(c)
            p5 = add(base5, ci(-5*c))
            live8 = count_translated_row(prepared, c, D)
            if (L+T+N+b+c) % 331 == 0:
                states_c = s4
                for k in range(5,9):
                    states_c = row_step(states_c, poly, k, D)
                assert live8 == len(states_c), (core,b,c)
                counts["translation_crosschecks"] += 1
            counts["rows8"] += live8
            if live8:
                counts["row8_core_bc"] += 1
                record("row8_b", b)
                record("row8_c", c)
                if final:
                    powers5 = [*powers,p4,p5]
                    tail = tail_box(poly, powers5, D)
                    if (L+T+N+b+c) % 6007 == 0:
                        pairs = []
                        for d in range(-41,42):
                            poly[6] = ci(d)
                            s6 = newton_next(poly, powers5, 6, D)
                            if not strict_house(s6, *BOUND[6], D):
                                continue
                            for e in range(-9,10):
                                poly[7] = ci(e)
                                s7 = newton_next(poly, [*powers5,s6], 7, D)
                                s8 = newton_next(poly, [*powers5,s6,s7], 8, D)
                                if strict_house(s7,*BOUND[7],D) and strict_house(s8,*BOUND[8],D):
                                    pairs.append((d,e))
                        poly[6] = poly[7] = ZERO
                        expected = ((len(pairs),min(d for d,e in pairs),max(d for d,e in pairs),
                                     min(e for d,e in pairs),max(e for d,e in pairs)) if pairs else None)
                        assert expected == tail, (core,b,c,expected,tail)
                        counts["tail_crosschecks"] += 1
                    if tail is not None:
                        n, dlo, dhi, elo, ehi = tail
                        counts["full_mothers"] += n
                        counts["active_tail_groups"] += 1
                        for key, value in (("b",b),("c",c),("d",dlo),("d",dhi),("e",elo),("e",ehi)):
                            record("full_"+key, value)
                        sha.update(json.dumps((b,c,*tail), separators=(",",":")).encode()+b"\n")
        poly[5] = ZERO
    return counts, ranges, sha.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max-row", type=int, choices=(7,8), default=8)
    ap.add_argument("--final", action="store_true")
    a = ap.parse_args()
    assert not a.out.exists()
    assert not a.final or a.max_row == 8
    entries = cores()
    totals, extrema, payload = Counter(), {}, hashlib.sha256()
    t_live_by_l = Counter()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        jobs = ((core,a.max_row,a.final) for core in entries)
        for core, (counts, ranges, h) in zip(entries, pool.map(audit_core, jobs, chunksize=16), strict=True):
            totals.update(counts)
            if counts["row7_core_b"]:
                t_live_by_l[core[0]] += 1
            for key,(x,y) in ranges.items():
                current = extrema.setdefault(key,[x,y])
                current[:] = min(current[0],x), max(current[1],y)
            payload.update(json.dumps((*core,h), separators=(",",":")).encode()+b"\n")
    summary = dict(status="complete", max_row=a.max_row, full_newton_mothers=a.final,
                   cores=len(entries), counts=dict(totals), actual_ranges=extrema,
                   row7_live_cores_by_L=dict(sorted(t_live_by_l.items())),
                   core_digest=payload.hexdigest())
    a.out.write_text(json.dumps(summary,sort_keys=True,indent=2)+"\n")
    print(json.dumps(summary,sort_keys=True))


if __name__ == "__main__":
    main()
