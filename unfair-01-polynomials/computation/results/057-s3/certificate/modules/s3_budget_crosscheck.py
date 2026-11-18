#!/usr/bin/env python3
"""Arbitrary-precision fallback and rowwise audit of native S3 budget trees."""
import argparse
import hashlib
import json
from collections import Counter
from math import isqrt
from pathlib import Path


def allowed(base, v, D, k, reverse):
    if k == (3 if reverse else 5):
        return None if v else (-base//2,-base//2)
    if not reverse and k==2:
        return (0,0) if (base,v) in ((0,0),(2,0)) else None
    if not v:
        assert base%2==0
        return (-base//2,1-base//2)
    radical = isqrt(D*v*v)
    if v<0:
        radical=-radical-1
    q=-((base+radical)//2)
    return (q,q)


def arbitrary(row, reverse):
    L,M,T,N,D,b,c,d,e=row
    F=[1,L,M,0,b,c,d,e,1]
    pos=5 if reverse else 3
    if reverse:F.reverse()
    H=[0]*17
    for j,x in enumerate(F):
        for k,y in enumerate(F):H[j+k]+=x*y
        H[j+pos]+=T*x
    H[2*pos]+=N
    p=[2*x for x in F]
    p[pos]=T
    states=[((1,),(1,),(1,))]
    max_live=0
    for k in range(1,97):
        next_states=[]
        for G,Q,W in states:
            U=sum(p[j]*G[k-j] for j in range(1,min(8,k)+1))
            V=sum(p[j]*Q[k-j] for j in range(1,min(8,k)+1))
            gg=allowed(U,G[k-pos] if k>=pos else 0,D,k,reverse)
            qq=allowed(V,-Q[k-pos] if k>=pos else 0,D,k,reverse)
            if gg is None or qq is None:
                continue
            internal=sum(x*y for x,y in zip(G[1:],Q[:0:-1],strict=True))
            carry=sum(H[j]*W[k-j] for j in range(1,min(16,k)+1))
            for target in (0,1):
                GQ=target-carry-internal
                lo,hi=max(gg[0],GQ-qq[1]),min(gg[1],GQ-qq[0])
                for g in range(lo,hi+1):
                    next_states.append((G+(g,),Q+(GQ-g,),W+(target-carry,)))
        if not next_states:
            return 'empty',k,max_live
        max_live=max(max_live,len(next_states))
        if max_live>10000:
            return 'state_cap',k,max_live
        states=next_states
    return 'depth_cap',96,max_live


def full_check(row):
    f,k,live=arbitrary(row,False)
    if f=='empty':
        return 'empty_forward',k,k,live
    r,h,lr=arbitrary(row,True)
    if r=='empty':
        return 'empty_reverse',k+h,h,lr
    return 'unresolved',k+h,0,0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--expected',type=int,default=869396)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    counts=Counter()
    hashes=[hashlib.sha256(),hashlib.sha256()]
    fallback=[]
    layers=0
    reverse_row=None
    with args.native.open('rb') as native, args.reference.open('rb') as reference:
        for original,checked in zip(native,reference,strict=True):
            hashes[0].update(original);hashes[1].update(checked)
            row,status,depth_f,depth_r,live,*_=json.loads(b'['+original.rstrip(b'\n')+b']')
            expected_row,wanted,expected_layers,death,expected_max,_=json.loads(checked)
            assert row==expected_row
            if status=='F':
                got=('empty_forward',depth_f,depth_f,live)
            else:
                got=full_check(row)
                fallback.append({'mother':row,'native_status':status,'native_forward':depth_f,
                                 'native_reverse':depth_r,'arbitrary':got})
            assert got == (wanted,expected_layers,death,expected_max),(row,got,(wanted,expected_layers,death,expected_max))
            counts[got[0]]+=1
            layers+=got[1]
            if got[0]=='empty_reverse':
                assert reverse_row is None
                reverse_row=row
    assert sum(counts.values())==args.expected and not counts.get('unresolved')
    result={'status':'complete','mothers':args.expected,'states':dict(counts),'layers':layers,
            'native_sha256':hashes[0].hexdigest(),'reference_sha256':hashes[1].hexdigest(),
            'fallbacks':fallback, 'reverse_only':reverse_row}
    with args.out.open('x') as f:
        f.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':
    if not __debug__:raise SystemExit('run without -O')
    main()
