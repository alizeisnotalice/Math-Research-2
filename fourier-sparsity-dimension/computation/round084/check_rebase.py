#!/usr/bin/env python3
"""Integer check for the sharp four-point rebase bound and its actual lifts."""
from collections import Counter
from itertools import combinations
import json
import platform
import random
import time

import phase_lifts as phase


def census(n):
    histogram=Counter(); winner=None
    for a,b in combinations(range(1,1<<n),2):
        c=a^b
        if c<b:
            continue
        flat=(0,a,b,c)
        for x in range(1,1<<n):
            P=tuple(sorted(u^x for u in flat))
            if P[0]!=x:
                continue
            total=sum(4**phase.pop(u) for u in P)
            histogram[total]+=1
            if total==40 and winner is None:
                winner=P
    theoretical=(((1<<n)-1)*((1<<n)-2)//6)*((1<<(n-2))-1)
    if sum(histogram.values())!=theoretical or min(histogram)!=40:
        raise AssertionError('unexpected nonzero affine plane minimum')
    return dict(n=n,affine_planes=theoretical,minimum_sum=40,minimizers=histogram[40],first_minimizer=winner)


def case(name):
    truth,n=phase.fixture_truth(name)
    old=phase.normalized_seed(truth,n)
    open_rows=list(phase.greedy_select(old,kmax=1)[0]['rows'])
    first,dense=phase.lift(old,open_rows,'disjoint')
    A=open_rows[0]
    rebased=[(A[0]^(1<<n),A[1]),(1<<n,0),(1<<(n+1),0)]
    second,dense2=phase.lift(first,rebased,'disjoint')
    direct=[phase.state_summary(s) for s in (old,first,second)]
    m0=phase.greedy_select(old,kmax=1)[0]['active']
    if len(second['coeff'])!=len(old['coeff'])-m0+40*m0 or len(first['coeff'])!=len(old['coeff'])+3*m0:
        raise AssertionError('sharp rebase does not multiply current tail by ten')
    return dict(seed=name,stages=direct,dense_checked=[True,dense,dense2],rows=[open_rows,rebased],
                quotient_jump=[1,3,9],activation_multiplier=10)


def random_linear_basis(rows,rng):
    new=[]; piv=[]; k=len(rows)
    while len(new)<k:
        v=rng.getrandbits(k)
        if phase.rank([*piv,v])==len(piv):
            continue
        piv.append(v)
        a=b=0
        for i,(p,q) in enumerate(rows):
            if v>>i&1:
                a^=p;b^=q
        new.append((a,b))
    return new


def random_rebase_checks():
    rng=random.Random(840831)
    counts=Counter()
    for name in ('rare-2','rare-3','rare-4','rare-5','r49','address-two'):
        seed=phase.normalized_seed(*phase.fixture_truth(name))
        opening=phase.greedy_select(seed,kmax=2)[0]
        for _ in range(4):
            old=seed; rows=opening['rows']; base=len(seed['coeff'])-opening['active']
            for i in range(4):
                rows=random_linear_basis(rows,rng)
                mode=('common','disjoint')[rng.randrange(2)]
                total=phase.phase_data(old,rows,mode)[-1]
                if total>4096:
                    break
                before=len(old['coeff'])-base
                v=old['n']
                old,dense=phase.lift(old,rows,mode)
                ratio=4 if mode=='common' or not i else 10
                if len(old['coeff'])-base<ratio*before:
                    raise AssertionError('a rebased step violated the tail growth')
                rows.extend((1<<j,0) for j in range(v,old['n']))
                counts[(mode,i,dense)]+=1
    return [dict(mode=m,stage=i+1,dense_checked=d,count=c) for (m,i,d),c in sorted(counts.items())]


if __name__=='__main__':
    t=time.perf_counter()
    data=dict(status='pass',host=platform.node(),python=platform.python_version(),
              affine_plane_census=[census(i) for i in range(3,8)],
              sharp_examples=[case(n) for n in ('rare-4','rare-5','nested-r49')],
              random_rebase_checks=random_rebase_checks())
    data['wall_seconds']=time.perf_counter()-t
    print(json.dumps(data,indent=2,sort_keys=True))
