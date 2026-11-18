#!/usr/bin/env python3
"""Recompute actual physical child truth tables for the strongest matching examples."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
import platform
import time

import matching_queries as match


def real_children(case,rows):
    n=case['n']
    pivots=[1<<(p.bit_length()-1) for p in match.basis(rows,n)]
    pivmask=sum(pivots)
    reps=[x for x in range(1<<n) if not x&pivmask]
    translated=match.span(pivots)
    truth=list(match.actual_truth(case))
    direct=[]
    for answer in range(1<<len(rows)):
        values=[]
        for z in reps:
            choices=[z|p for p in translated if all(match.parity((z|p)&q)==((answer>>i)&1)
                                                     for i,q in enumerate(rows))]
            if len(choices)!=1:
                raise AssertionError(('invalid physical coordinate chart',rows,answer,z))
            values.append(truth[choices[0]])
        W=match.fwht(values)
        direct.append({u:w for u,w in zip(reps,W) if w})
    return direct


def verify_real(case,rows):
    direct=real_children(case,rows)
    exact=match.restrictions(case,rows)
    normalizer=1<<(case['n']-len(rows))
    for f,g in zip(direct,exact):
        if set(f)!=set(g) or any(case['denom']*f[q]!=normalizer*g[q] for q in f):
            raise AssertionError((case['name'],rows,'actual direct truth transform disagrees'))
    return dict(H_basis=list(rows),local_type=list(match.queries_type(rows,case['address_dim'],case['n'])),
                child_pairs=match.child_pairs(case,direct),
                child_integer_walsh=direct,
                child_P_squared=[list((z.numerator,z.denominator)) for p in match.child_pairs(case,direct) for z in [match.pair_score(p)]])


def through_line(case,q):
    ds=set(match.directions(case))
    hist=Counter(); witnesses={}
    match.one_diff_planes(case,q,ds,match.restrictions(case,(q,)),hist,witnesses)
    for r in sorted({min(w,w^q) for w in ds if w!=q}):
        rows=(q,r)
        m=2+(r in ds and r^q in ds)
        a,b=match.queries_type(rows,case['address_dim'],case['n'])
        match.summarize_type(hist,witnesses,a,b,m,match.child_pairs(case,match.restrictions(case,rows)),rows)
    if sum(hist.values())!=(1<<(case['n']-1))-1:
        raise AssertionError('not all planes through the prescribed line were checked')
    overall=max(hist,key=lambda k:max(map(match.pair_score,k[3])))
    parent=match.pair_score((len(case['coeff']),case['R']))
    recovery=sum(c for k,c in hist.items() if max(map(match.pair_score,k[3]))>=parent)
    strict=sum(c for k,c in hist.items() if max(map(match.pair_score,k[3]))>parent)
    return dict(fixed_direction=q,n_planes=sum(hist.values()),n_recovering=recovery,n_improving=strict,
                exact_classes=[dict(address_intersection=a,data_intersection=b,active_differences=m,child_pairs=p,
                                    n_kernels=c) for (a,b,m,p),c in sorted(hist.items())],
                best=verify_real(case,witnesses[overall]))


def pure_address_planes():
    return sorted({tuple(sorted((x,y,x^y))[:2]) for x,y in combinations(range(1,8),2)})


def main():
    begin=time.perf_counter()
    refs={z['name']:z for f in ('matching-mini.json','matching-macbook.json','guarded-odd-mini.json')
          for z in json.load(open(f))['cases']}
    cases=[]
    for name,p in refs.items():
        case=match.build_fixture(name)
        checks={}
        for level in p['levels']:
            for group in [level['overall_best'],*[x['best'] for x in level['types']]]:
                if group['H_basis'] is not None:
                    h=tuple(group['H_basis'])
                    if h in checks:
                        continue
                    # Direct transforms for the absolute optima and the guarded R49 near ties.
                    if group!=level['overall_best'] and (name!='r49:private-guarded-bent' or group['p2'][0]*13<49*group['p2'][1]):
                        continue
                    checks[h]=verify_real(case,h)
        extra=None
        if name in ('r49:private-guarded-bent','r49:private-guarded-odd'):
            for h in pure_address_planes():
                if h not in checks:
                    checks[h]=verify_real(case,h)
                if max(map(match.pair_score,checks[h]['child_pairs']))!=Fraction(*p['P_squared']):
                    raise AssertionError('pure two-address witness unexpectedly differs')
        if name=='r49:private-guarded-odd':
            extra=through_line(case,66)
        cases.append(dict(name=name,real_witnesses=list(checks.values()),all_planes_through_q66=extra))
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),cases=cases,
                          wall_seconds=time.perf_counter()-begin),sort_keys=True))


if __name__=='__main__':
    main()
