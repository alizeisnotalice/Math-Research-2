#!/usr/bin/env python3
"""Locate and verify skew-rank gains when the sparse phase support has holes."""
from fractions import Fraction
import json
import platform
import random
import time

import phase_lifts as phase


def scan():
    seeds=[('cubic',5,4096),('rare',6,4096),('rare',7,1024)]
    best={}; counts={}
    for family,n,total in seeds:
        rng=random.Random(840830+100*n)
        for i in range(total):
            truth=phase.random_truth(n,family,rng)
            fg=phase.normalized_seed(truth,n)
            base=phase.state_summary(fg)
            if base['s']==1<<base['v']:
                continue
            pbase=Fraction(*base['P_squared'])
            for t in phase.greedy_select(fg):
                k=t['k']
                maxrank=max(map(int,t['disjoint_weight_histogram']))*2
                for key,flag,pfield in [('k3-disjoint',k>=3,'disjoint_p2'),('rank4-disjoint',maxrank>=4,'disjoint_p2'),
                                        ('rank6-disjoint',maxrank>=6,'disjoint_p2'),('k3-common',k>=3,'common_p2')]:
                    if not flag:
                        continue
                    value=Fraction(*t[pfield]); gain=value/pbase
                    if gain<=1:
                        continue
                    counts[key]=counts.get(key,0)+1
                    if key not in best or gain>Fraction(*best[key]['gain_squared']):
                        best[key]=dict(case=key,family=family,truth=truth,n=n,choice=t,p2_before=base['P_squared'],
                                       p2_after=t[pfield],gain_squared=[gain.numerator,gain.denominator],max_polar_rank=maxrank)
    for key,record in best.items():
        s=phase.normalized_seed(record['truth'],record['n'])
        mode='common' if key.endswith('common') else 'disjoint'
        record['explicit_trace']=phase.follow_recursive(s,record['choice'],mode)
        record['actual_first_lift']=phase.state_summary(phase.lift(s,record['choice']['rows'],mode)[0])
    return dict(inputs=[dict(kind=k,n=n,tables=c) for k,n,c in seeds],selected_gains=counts,witnesses=list(best.values()))


def arbitrary_rank_fixture(k=3,core=7):
    # On an affine address half-flat use a full-spectrum pointflip. On the other half,
    # four singleton addresses select the first three private characters and z1+z2.
    # Use actual Boolean truth to validate the non-full-spectrum spectrum.
    t=3; n=t+core+k
    leaves=(1,2,4,3)
    values=[]
    for x in range(1<<n):
        a=x&7
        if a<4:
            flip=phase.parity((x>>(t+core))&leaves[a])
        else:
            flip=int((x>>t)&((1<<core)-1)==0)
        values.append(1-2*flip)
    W=phase.fwht(values)
    coeff={u:v for u,v in enumerate(W) if v}
    state=dict(n=n,denom=1<<n,coeff=coeff)
    if phase.rank((1<<n)|u for u in coeff)!=n+1:
        raise AssertionError('poor seed affine coordinates')
    rows=[(1<<(t+core+i),0) for i in range(k)]
    mode='disjoint'
    child,checked=phase.lift(state,rows,mode,dense_cap=19)
    a=phase.state_summary(state); b=phase.state_summary(child)
    if Fraction(*b['P_squared'])<=Fraction(*a['P_squared']):
        raise AssertionError('chosen skew-rank example is not an improvement')
    return dict(k=k,core=core,seed=a,child=b,dense_child_transform_checked=checked,rows=rows,
                label_weight_histogram=dict(__import__('collections').Counter(phase.pop(phase.evaluate_affines(rows,u)) for u in coeff)))


if __name__=='__main__':
    start=time.perf_counter()
    data=scan()
    data['dense_fallback_rank4']=arbitrary_rank_fixture()
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),
                          total_wall_seconds=time.perf_counter()-start,**data),indent=2,sort_keys=True))
