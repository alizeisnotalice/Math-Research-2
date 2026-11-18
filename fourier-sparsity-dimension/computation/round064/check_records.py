"""Compare a function to all affine restrictions before testing the record inequalities."""
from fractions import Fraction
from itertools import combinations
from collections import Counter
import json
import platform


def independent(points):
    basis=[]
    span={0}
    for p in points:
        if p not in span:
            basis.append(p)
            span |= {x^p for x in span}
    return basis


def points(basis):
    output=[0]
    for b in basis:
        output += [x^b for x in output]
    return output


def spaces(n):
    found={frozenset({0})}
    for j in range(1,n+1):
        for basis in combinations(range(1,1<<n),j):
            found.add(frozenset(points(basis)))
    return sorted(found,key=lambda s:(len(s),tuple(sorted(s))))


def spectrum(values):
    row=values[:]
    step=1
    while step<len(row):
        for k in range(0,len(row),step*2):
            for j in range(k,k+step):
                a,b=row[j],row[j+step]
                row[j],row[j+step]=a+b,a-b
        step*=2
    return [i for i,v in enumerate(row) if v]


def check(n,masks):
    subspaces=spaces(n)
    affine=[]
    for V in subspaces:
        ordered=points(independent(sorted(V)))
        for x in range(1<<n):
            if x==min(x^v for v in V):
                affine.append([x^v for v in ordered])
    checked=records=0
    for mask in masks:
        f=[1-2*((mask>>x)&1) for x in range(1<<n)]
        S=spectrum(f)
        d=len(independent(S))
        if not d:
            continue
        if any(len(independent(T))**2*len(S)>d*d*len(T)
               for A in affine if (T:=spectrum([f[x] for x in A]))):
            continue
        records+=1
        for H in subspaces:
            B=independent(sorted(H))
            ambient=independent(B+list(range(1<<n)))
            values=points(ambient)
            quotient={x:i>>len(B) for i,x in enumerate(values)}
            counts=Counter(quotient[x] for x in S)
            I=independent(sorted(counts,key=lambda k:counts[k]))
            for j in range(n+2):
                M=1<<j
                cap=sum(min(k,M) for k in counts.values())
                lower=max(0,len(I)-sum(Fraction(counts[z]-1,2*M) for z in I))
                assert Fraction(cap,len(S))>=Fraction(lower*lower,d*d)
                checked+=1
    return dict(n=n,tested=len(masks),affine_domains=len(affine),records=records,
                simultaneous_direction_scale_checks=checked)


print(json.dumps(dict(python=platform.python_version(),results=[
    check(3,range(256)),check(4,[i*257 for i in range(256)]+[i*241 for i in range(256)])]),indent=2))
