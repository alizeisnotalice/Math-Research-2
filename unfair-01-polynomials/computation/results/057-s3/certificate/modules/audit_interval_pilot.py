#!/usr/bin/env python3
"""Verify S3 pilot survivors and endpoint-aware deterministic sample with exact norm Schur."""
from pathlib import Path
import argparse
from collections import Counter
from math import gcd
from itertools import chain
import hashlib
import json

POINTS=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))


def norm(row):
    L,M,T,N,D,b,c,d,e=row
    F=[1,L,M,0,b,c,d,e,1]
    result=[0]*17
    for j,a in enumerate(F):
        for k,b in enumerate(F):
            result[j+k]+=a*b
        result[j+3]+=T*a
    result[6]+=N
    return result


def strict_disk(p):
    while len(p)>1:
        lead, tail=p[0], p[-1]
        if abs(lead)<=abs(tail):
            return False
        p=[lead*x-tail*y for x,y in zip(p[:-1],p[:0:-1],strict=True)]
        scale=gcd(*p)
        p=[x//scale for x in p]
    return True


def exact_status(row):
    L,M,T,N,D,b,c,d,e=row
    F=[1,L,M,0,b,c,d,e,1]
    for u,v in POINTS:
        f=sum(c*u**k*v**(8-k) for k,c in enumerate(F))
        irr=u**3*v**5
        lhs=2*f+T*irr
        if lhs<=0 or lhs*lhs<=D*irr*irr:
            return 'positive_reject'
    h=norm(row)
    if not strict_disk([a*2**k for k,a in enumerate(h)][::-1]):
        return 'outer_reject'
    if not strict_disk([a*2**(16-k) for k,a in enumerate(h)]):
        return 'inner_reject'
    return 'pass'


def selected(group):
    L,T,N,D,b,c,dlo,dhi,elo,ehi,zlo,zhi,num=group
    chosen=[dlo,(dlo+dhi)//2,dhi]
    for d in chosen:
        a=max(elo,zlo+L*d)
        e=min(ehi,zhi+L*d)
        if a>e:
            continue
        for val in (a,(a+e)//2,e):
            yield L,L*L,T,N,D,b,c,d,val


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=Path,required=True)
    parser.add_argument('--pass-file',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    published={tuple(json.loads(line)) for line in args.pass_file.read_bytes().splitlines()}
    stats=Counter()
    for row in published:
        assert exact_status(row)=='pass', row
    touched=set()
    for line in args.groups.read_text().splitlines():
        for row in selected(tuple(map(int,line.split()))):
            if row in touched:continue
            touched.add(row)
            result=exact_status(row)
            stats[result]+=1
            assert (result=='pass')==(row in published), (row,result)
    data={'published_pass_checked':len(published),'sample_unique':len(touched),
          'sample_counts':dict(stats),'arithmetic':'Python integer norm and exact integer Schur',
          'pass_file_sha256':hashlib.sha256(args.pass_file.read_bytes()).hexdigest(),
          'groups_sha256':hashlib.sha256(args.groups.read_bytes()).hexdigest(),'status':'complete'}
    with args.out.open('x') as stream:
        stream.write(json.dumps(data,sort_keys=True,indent=2)+'\n')
    print(json.dumps(data,sort_keys=True))


if __name__=='__main__':
    main()
