#!/usr/bin/env python3
"""Exact algebraic rectangle-cell checks; no solver, random search or outside halo."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import time
import sympy as sp


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    assert __debug__, 'assertions required'
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    start = time.perf_counter()
    inp = json.loads(Path(args.input).read_text())
    half = sp.Rational(1, 2)
    quarter = sp.Rational(1, 4)
    r1 = sp.Rational(3, 4)
    r2 = 1-sp.sqrt(2)/4
    rows = []
    for text in inp['exact_geometry_fixtures']:
        t = sp.sympify(text)
        pts = [(0,0),(t,0),(t,t),(0,t)]
        def e(mask, radius):
            pp = [pts[i] for i in range(4) if mask >> i & 1]
            return tuple(map(sp.sympify, (max(p[0] for p in pp)-radius,
                min(p[0] for p in pp)+radius,
                max(p[1] for p in pp)-radius,min(p[1] for p in pp)+radius)))
        def rho(k):
            return r1 if k == 1 else 1-sp.sqrt(sp.Rational(k,4))/2
        internal = [e(1 << i, quarter) for i in range(4)]+[e(15,half)]
        ers = {m:e(m,half) for m in range(1,15)}
        ks = {k:(t-rho(k),rho(k),t-rho(k),rho(k)) for k in (1,2,3)}
        c = (-half,t+half,-half,t+half)
        allrect = internal+list(ers.values())+list(ks.values())+[c]
        xs = sorted(set(z for r in allrect for z in r[:2]))
        ys = sorted(set(z for r in allrect for z in r[2:]))
        xi, yi = {x:i for i,x in enumerate(xs)}, {y:i for i,y in enumerate(ys)}
        def cells(r):
            a,b,c_,d=r
            if a >= b or c_ >= d:
                return set()
            return set(itertools.product(range(xi[a],xi[b]),range(yi[c_],yi[d])))
        def area(cc):
            return sp.simplify(sum((xs[i+1]-xs[i])*(ys[j+1]-ys[j]) for i,j in cc))
        I = set().union(*(cells(r) for r in internal))
        carrier = {m:cells(ers[m])-cells(ks[m.bit_count()]) for m in ers}
        singleton = set().union(*(carrier[1<<i] for i in range(4)))
        pair_masks = [m for m in ers if m.bit_count()==2 and m not in (5,10)]
        pair = set().union(*(carrier[m] for m in pair_masks))
        k1 = cells(ks[1])
        inside_pair = pair & k1
        pure_eligible=[]
        for m in range(1,16):
            pp=[pts[i] for i in range(4) if m>>i&1]
            d=max(max(p[j] for p in pp)-min(p[j] for p in pp) for j in (0,1))
            if sp.Rational(m.bit_count(),4)>d*d:
                pure_eligible.append(m)
        assert pure_eligible == [1,2,4,8,15]
        assert singleton == cells(c)-k1
        assert sp.simplify(area(I & k1)-5*(1-t)**2)==0
        assert sp.simplify(area(inside_pair)-(1-t)*(sp.sqrt(2)-1))==0
        assert not inside_pair & I
        strips = [carrier[m]&k1 for m in pair_masks]
        assert all(not x&y for x,y in itertools.combinations(strips,2))
        for m in ers:
            if m.bit_count()==3 or m in (5,10):
                assert cells(ers[m]) <= I
        total = I | set().union(*carrier.values())
        F = sp.Rational(15,4)-5*t+5*t*t+(sp.sqrt(2)-1)*(1-t)
        assert sp.simplify(area(total)-F)==0
        rows.append({'t':str(t),'pure_eligible_masks':pure_eligible,
            'proper_captured_sets_checked':14,'cells':len(xs[:-1])*len(ys[:-1]),
            'internal_halo_area':str(area(I)),
            'singleton_carrier_area':str(area(singleton)),
            'I_intersect_K1_area':str(area(I&k1)),
            'new_adjacent_pair_area':str(area(inside_pair)),
            'coarse_union_area':str(area(total)),
            'coarse_union_decimal_informational':str(sp.N(F,20)), 'status':'PASS'})
    t=sp.symbols('t', real=True)
    F=sp.Rational(15,4)-5*t+5*t*t+(sp.sqrt(2)-1)*(1-t)
    assert sp.simplify(sp.diff(F,t)-(10*t-4-sp.sqrt(2)))==0
    # derivative >= 5sqrt3-4-sqrt2 > 0, using sqrt3>17/10,sqrt2<3/2.
    assert sp.Rational(17,10)**2<3 and sp.Rational(3,2)**2>2
    assert 5*sp.Rational(17,10)-4-sp.Rational(3,2)>0
    top=sp.simplify(F.subs(t,sp.Rational(9,10)))
    assert sp.simplify(top-(32+sp.sqrt(2))/10)==0
    assert sp.Rational(3,2)**2>2  # top < 67/20 < 7/2
    result={'status':'PASS','scope':inp['purpose'], 'input_sha256':digest(args.input),
        'source_sha256':digest(__file__),'python':platform.python_version(),'sympy':sp.__version__,
        'pid':os.getpid(),'fixtures':rows,'all_t_formula':str(F),
        'uniform_upper':str(top),'strict_rational_upper':'67/20',
        'strict_threshold':'7/2','elapsed_seconds':time.perf_counter()-start}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','fixtures':len(rows),'upper':str(top),
                      'elapsed_seconds':result['elapsed_seconds']}))

if __name__=='__main__':
    main()
