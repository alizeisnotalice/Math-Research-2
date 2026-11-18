"""Exact Walsh spectra and truth checks for nested tail copies (Python 3.10+)."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
import platform
import random


def rank(points):
    pivots={}
    for a in points:
        while a and a.bit_length() in pivots:
            a ^= pivots[a.bit_length()]
        if a:
            pivots[a.bit_length()]=a
    return len(pivots)


def walsh(table):
    row=table.copy()
    w=1
    while w<len(row):
        for i in range(0,len(row),2*w):
            for j in range(i,i+w):
                a,b=row[j],row[j+w]
                row[j],row[j+w]=a+b,a-b
        w*=2
    return {i:c for i,c in enumerate(row) if c}


def lift(state,h,B):
    n,den,spec=state
    mask=(1<<h)-1
    address=h+B*(n-h)
    new={}
    for x,c in spec.items():
        a,b=x&mask,x>>h
        if not b:
            new[a]=c*B
            continue
        for leaf in range(B):
            node=B+leaf
            terms=[(0,1)]
            while node>1:
                bit=1<<(address+(node//2)-1)
                sign=1 if node%2==0 else -1
                terms += [(z|bit,sign*v) for z,v in terms]
                node//=2
            origin=a|(b<<(h+leaf*(n-h)))
            for z,v in terms:
                assert origin|z not in new
                new[origin|z]=v*c
    r=rank(x>>h for x in spec)
    assert r
    q=sum(x>>h!=0 for x in spec)
    assert len(new)==len(spec)+(B**2-1)*q
    assert rank(new)==rank(spec)+(B-1)*(r+1)
    return address+B-1,den*B,new


def unlift(z,n,h,B):
    address=h+B*(n-h)
    node=1
    while node<B:
        node=2*node+((z>>(address+node-1))&1)
    return (z&((1<<h)-1))|(((z>>(h+(node-B)*(n-h)))&((1<<(n-h))-1))<<h)


def spectrum_eval(spec,z):
    return sum(-v if (i&z).bit_count()%2 else v for i,v in spec.items())


def hull_matroid_tests():
    # Compute principal flats and verify the fractional base inequalities
    # directly on every flat in the effective span.
    subs={frozenset([0])}
    for i in range(1,16):
        subs |= {s|frozenset(x^i for x in s) for s in subs}
    spectra={frozenset(walsh([1-2*((t>>i)&1) for i in range(16)]))-{0}
             for t in range(65536)}
    biggest=Fraction(0)
    histogram=Counter()
    for E in spectra:
        h=frozenset([0])
        x={}
        flats=[s for s in subs if rank(s&E)==rank(s)]
        levels=0
        while not E<=h:
            cand=[s for s in flats if h<s]
            def slope(s):
                return Fraction(rank(s)-rank(h),len((s-h)&E))
            a=min(slope(s) for s in cand)
            chosen=max((s for s in cand if slope(s)==a),key=len)
            for b in (chosen-h)&E:
                x[b]=a
            h=chosen
            levels+=1
        assert sum(x.values())==rank(E)
        assert all(sum((x[b] for b in E&s),Fraction(0))<=rank(s) for s in subs)
        cost=sum((z*z for z in x.values()),Fraction(0))
        biggest=max(biggest,cost)
        histogram[levels]+=1
    return dict(distinct_nonzero_supports=len(spectra),maximum_squared_cost=str(biggest),
                principal_level_counts=dict(sorted(histogram.items())))


def run():
    rng=random.Random(7401)
    checks=full=0
    for t in range(256):
        row=[1-2*((t>>i)&1) for i in range(8)]
        spec=walsh(row)
        for h in (0,1,2):
            if not any(x>>h for x in spec):
                continue
            for B in (2,4):
                state=lift((3,8,spec),h,B)
                n,den,new=state
                values=range(1<<n) if n<=7 else [rng.getrandbits(n) for _ in range(12)]
                for z in values:
                    assert spectrum_eval(new,z)==den*row[unlift(z,3,h,B)]
                if n<=7:
                    full+=1
                checks+=1
    # Full-support three-variable sign flip. Weights by shells: 1,2,4.
    row=[-1]+[1]*7
    state=lift((3,8,walsh(row)),2,2)
    n=state[0]
    final_n,den,new=lift(state,1,2)
    assert len(new)==1+1+4*2+16*4 and rank(new)==1+2+4+4-1
    for _ in range(500):
        z=rng.getrandbits(final_n)
        assert spectrum_eval(new,z)==den*row[unlift(unlift(z,n,1,2),3,2,2)]
    return dict(python=platform.python_version(),seed=7401,one_stage=checks,
                exhaustive_physical_checks=full,two_stage_s=len(new),
                two_stage_d=rank(new),matroid=hull_matroid_tests())


if __name__=='__main__':
    print(json.dumps(run(),indent=2))
