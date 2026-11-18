"""Exact pointwise denominator and relative rank checks, all four-bit tables."""
from collections import Counter
import json

def dimension(columns):
    basis={}
    for v in columns:
        while v and v.bit_length() in basis:
            v^=basis[v.bit_length()]
        if v: basis[v.bit_length()]=v
    return len(basis)

def main():
    # Group by support, keeping the largest observed exponent at each frequency.
    masks=[sum(1<<x for x in range(16) if (x&a).bit_count()%2) for a in range(16)]
    supports={}
    for truth in range(65536):
        w=[16-2*(truth^mask).bit_count() for mask in masks]
        s=tuple(i for i,z in enumerate(w) if z)
        exponents=supports.setdefault(s,[0]*16)
        for i in s:
            lowbit=w[i]&-w[i]
            exponents[i]=max(exponents[i],5-lowbit.bit_length())
    spaces={frozenset({0})}
    for _ in range(4):
        spaces |= {space|{x^z for x in space} for space in spaces for z in range(16)}
    flats={frozenset(x^z for x in space) for space in spaces for z in range(16)}
    results=Counter()
    for s,exponent in supports.items():
        for flat in flats:
            t=5-len(flat).bit_length()
            xs=set(s)&flat
            nullity=len(xs)-dimension(x+16 for x in xs)
            for x in xs:
                k=max(0,exponent[x]+1-t)
                assert nullity>=(1<<k)-1-k
                results[t]+=1
    return dict(tables=65536,supports=len(supports),flats=len(flats),label_checks=results)
print(json.dumps(main(),indent=2))
