#!/usr/bin/env python3
"""Exact test of target-address composition, including balanced seeds."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from exact_baseline import walsh,basis


def check(v,a,truth_check=True):
    n=(len(v)-1).bit_length();t=1 << a;w=walsh(v);N=len(v)
    support={}
    if w[0]:
        support[0]=t*w[0]
    for i in range(t):
        for u in range(1,N):
            if not w[u]:
                continue
            for z in range(t):
                key=(u << (a+i*n))|z
                assert key not in support
                support[key]=w[u]*((-1)**((i&z).bit_count()%2))
    s0=sum(x!=0 for x in w);d0=len(basis(i for i,x in enumerate(w) if x));delta=int(w[0]!=0)
    assert len(support)==delta+t*t*(s0-delta)
    assert len(basis(support))==t*d0+a
    if truth_check:
        truth=[v[(x >> (a+n*(x&(t-1))))&(N-1)] for x in range(1 << (a+n*t))]
        direct=walsh(truth)
        multiplier=len(truth)//(N*t)
        assert {i:z for i,z in enumerate(direct) if z}=={i:z*multiplier for i,z in support.items()}
    return {'n':n,'a':a,'seed_s':s0,'seed_d':d0,'delta_zero':delta,
            's':len(support),'d':len(basis(support)),'truth_check':truth_check}


ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
counts=[];witnesses=[]
for n in range(1,4):
    N=1 << n;count=0
    for mask in range(1,(1 << N)-1):
        v=[1-2*((mask >> x)&1) for x in range(N)]
        r=check(v,1);count+=1
        if mask in (1,3,(1 << (N//2))-1):
            witnesses.append(check(v,2))
    counts.append({'seed_n':n,'all_nonconstant_seeds':count,'target_arity':2})
sequences=[]
for arities in ([2]*12,[4]*8,[2,4,8,16,2,8,4],[16,2,2,4,2]):
    P=1;weighted=Fraction(0)
    for t in arities:
        a=t.bit_length()-1;P*=t;weighted+=Fraction(a,P)
        assert weighted<=1-Fraction(1,P)
    sequences.append({'arities':arities,'product':P,'rank_increment_over_product':str(weighted),
                      'upper_bound':str(1-Fraction(1,P))})
out={'schema':'a28-target-address-v1','python':sys.version,'arithmetic':'integer and Fraction',
     'seed':None,'tolerance':0,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'counts':counts,'arity4_witnesses':witnesses,'variable_arity_checks':sequences}
args.out.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'passed','arity2_seeds':sum(x['all_nonconstant_seeds'] for x in counts),
                  'arity4_witnesses':len(witnesses)}))
