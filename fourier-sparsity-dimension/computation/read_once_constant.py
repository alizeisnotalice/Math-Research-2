#!/usr/bin/env python3
"""Exact structural DP for fresh-node selector trees with constant leaves.

States are (spectral sparsity, absolute mean, affine-full flag), indexed by essential selector
count. Global sign normalization and swapping disjoint child blocks preserve
these statistics. No floating-point arithmetic is used for proof checks.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
from exact_baseline import walsh, basis


def neg(tree):
    if isinstance(tree,int):
        return -tree
    return [neg(tree[0]),neg(tree[1])]


def table(tree):
    if isinstance(tree,int):
        return [tree]
    left,right=map(table,tree)
    # Root low bit, left and right child coordinates in separate blocks.
    return [left[x] if z==0 else right[y]
            for y in range(len(right)) for x in range(len(left)) for z in (0,1)]


def direct(v):
    return [sum(value*(-1 if (a & x).bit_count()%2 else 1)
                for x,value in enumerate(v)) for a in range(len(v))]


def audit_tree(tree,d,s,mu,epsilon):
    v=table(tree)
    assert len(v)==1 << d
    w=walsh(v)
    assert sum(a!=0 for a in w)==s
    assert len(basis(i for i,a in enumerate(w) if a))==d
    lifted=[i | (1 << d) for i,a in enumerate(w) if a]
    assert len(basis(lifted))==d+epsilon
    assert Fraction(sum(v),len(v))==mu
    assert walsh(w)==[len(v)*x for x in v]
    # Full exact g^2=1 via inverse-transform identity; direct characters at d<=8.
    direct_checked=d<=8
    if direct_checked:
        assert w==direct(v)
    return {'tree':tree,'truth_sha256':hashlib.sha256(bytes(x==1 for x in v)).hexdigest(),
            'walsh_sha256':hashlib.sha256(json.dumps(w).encode()).hexdigest(),
            'truth_entries':len(v),'direct_character_check':direct_checked}


def run(max_d):
    dp=[{(1,Fraction(1),1):1}]
    records=[]
    for d in range(1,max_d+1):
        states={}
        for dl in range((d-1)//2+1):
            dr=d-1-dl
            for (sl,ml,el),tl in dp[dl].items():
                for (sr,mr,er),tr in dp[dr].items():
                    for sign in (-1,1):
                        mright=sign*mr
                        mu=(ml+mright)/2
                        if dl==dr==0 and ml==mright:
                            continue  # constant root has no essential selector
                        s=2*(sl-int(bool(ml))+sr-int(bool(mr)))
                        s+=int(bool(ml+mright))+int(bool(ml-mright))
                        assert 4*s >= (d+1)**2+7*int(bool(mu))
                        epsilon=0 if dl==dr==0 else int(el or er)
                        assert 4*s >= (d+1+epsilon)**2
                        key=(s,abs(mu),epsilon)
                        if key not in states:
                            t=[tl,tr if sign==1 else neg(tr)]
                            states[key]=neg(t) if mu<0 else t
        dp.append(states)
        best=[]
        for balanced in (True,False):
            relevant=[k for k in states if (k[1]==0)==balanced]
            if relevant:
                s,mu,epsilon=min(relevant)
                best.append({'balanced':balanced,'s':s,'mu':str(mu),'epsilon':epsilon,
                             'exact_witness':audit_tree(states[s,mu,epsilon],d,s,mu,epsilon)})
        records.append({'d':d,'state_count':len(states),'minimum':best})
        print(json.dumps({'d':d,'states':len(states),'minima':[(r['balanced'],r['s']) for r in best]}),flush=True)
    return records


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--max-d',type=int,default=12)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();start=time.monotonic()
    records=run(a.max_d)
    data={'schema':'a28-read-once-constant-v1','max_d':a.max_d,'arithmetic':'integer and Fraction',
          'seed':None,'tolerance':0,'python':sys.version,'platform':platform.platform(),
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'elapsed_seconds':time.monotonic()-start,'records':records,
          'scope':'all structural (s,absolute mean,affine-full flag) states through max_d; exact truth tables only saved minima',
          'general_claim':'proved analytically in round003; finite DP is supporting verification'}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    main()
