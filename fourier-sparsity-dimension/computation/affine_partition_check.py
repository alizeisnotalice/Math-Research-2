#!/usr/bin/env python3
"""Exhaust all affine partitions of F2^3 and verify exact spectra and rank."""
import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
from exact_baseline import basis, walsh


def partitions(n):
    N=1 << n
    blocks=[]
    for mask in range(1,1 << N):
        xs=[i for i in range(N) if (mask >> i)&1]
        if len(xs)&(len(xs)-1):
            continue
        H={i ^ xs[0] for i in xs}
        if all(x ^ y in H for x in H for y in H):
            blocks.append((mask,xs,H))
    by_point=[[b for b in blocks if b[0] >> i & 1] for i in range(N)]
    def visit(left,chosen):
        if not left:
            yield chosen
            return
        p=(left & -left).bit_length()-1
        for b in by_point[p]:
            if b[0] & left == b[0]:
                yield from visit(left ^ b[0],chosen+[b])
    return blocks,visit((1 << N)-1,[])


def check(n):
    blocks,parts=partitions(n)
    histogram=Counter()
    checks=0
    repeated_checks=0
    for part in parts:
        t=len(part)
        normals=[]
        s=0
        for _,xs,H in part:
            W=[a for a in range(1 << n) if all((a & h).bit_count()%2==0 for h in H)]
            normals.extend(W)
            s+=len(W)
        R=len(basis(normals))
        assert R<=t-1
        assert sum(Fraction(len(xs),1 << n) for _,xs,H in part)==1
        assert s>=t*t
        for mode in ('independent','dependent','modulated-independent','modulated-dependent'):
            spec={}
            for i,(_,xs,H) in enumerate(part):
                b=(1 << i) if mode in ('independent','modulated-independent') else i
                mod=(3*i+1) % (1 << n) if mode.startswith('modulated') else 0
                v=[((-1)**((mod & x).bit_count()%2)) if x in xs else 0 for x in range(1 << n)]
                w=walsh(v)
                for a,c in enumerate(w):
                    if c:
                        freq=a | (b << n)
                        assert freq not in spec
                        spec[freq]=c
            d=len(basis(spec))
            assert len(spec)==s and d<=R+t and 4*s>=(d+1)**2
            if mode in ('independent','modulated-independent'):
                assert d==R+t
            checks+=1
        for all_same in (False,True):
            bands={}
            masses=Counter()
            multiplicities=Counter()
            for i,(_,xs,H) in enumerate(part):
                b=0 if all_same else i%2
                mod=(3*i+1) % (1 << n)
                bands.setdefault(b,[0]*(1 << n))
                for x in xs:
                    bands[b][x]=(-1)**(i%2+(mod & x).bit_count()%2)
                masses[b]+=len(xs)
                multiplicities[b]+=1
            support=[]
            for b,v in bands.items():
                w=walsh(v)
                supp=[i for i,c in enumerate(w) if c]
                assert len(supp)*masses[b]>=1 << n
                support.extend(a | (b << n) for a in supp)
            B=len(bands);K=max(multiplicities.values());d=len(basis(support))
            assert len(support)>=B*B and d<=2*t-1
            assert 4*K*K*len(support)>=(d+1)**2
            repeated_checks+=1
        histogram[t,R,s]+=1
    return {'n':n,'affine_blocks':len(blocks),'partitions':sum(histogram.values()),
            'spectral_checks':checks,'repeated_label_checks':repeated_checks,'histogram':[
                {'t':t,'R':R,'s':s,'count':c} for (t,R,s),c in sorted(histogram.items())]}


def ledger():
    out=[]
    for m in range(1,17):
        s=4**m;M=2**m;ell=s;cost=0;steps=0
        for j in range(m,0,-1):
            L=2**(m+j);k=4**j;delta=k//2;q=2**j;count=2**(m-j)
            assert ell==L
            for _ in range(count):
                assert k*s>=ell*ell and k<=ell
                assert 2*delta==k and q*q==k
                ell-=delta;cost+=q;steps+=1
            assert ell==L//2
        assert ell==M and cost==m*M
        out.append({'m':m,'s':s,'terminal_ell':ell,'cost':cost,'steps':steps})
    return out


ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
data={'arithmetic':'exact integers and rational numbers','random_seed':None,
      'partitions':[check(n) for n in range(4)],'scalar_ledger':ledger()}
a.out.write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'partitions':[{k:v for k,v in r.items() if k!='histogram'} for r in data['partitions']],
                  'ledger_max_m':16,'status':'passed'}))
