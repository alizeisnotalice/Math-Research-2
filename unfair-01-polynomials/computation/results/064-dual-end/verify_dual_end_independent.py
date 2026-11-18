#!/usr/bin/env python3
"""Independently check the finite arithmetic for the dual-end ramification bound."""

import argparse
from collections import Counter
from math import factorial, prod
from pathlib import Path
import json

import sympy as sp


def partitions(n, minimum=1):
    if n==0:
        yield ()
    for k in range(minimum,n+1):
        for left in partitions(n-k,k):
            yield k,*left


def shapes(cost):
    for parameters in partitions(cost):
        yield tuple(k+1 for k in parameters)


def cycle_profiles(sizes):
    if not sizes:
        yield ()
    else:
        for cluster in partitions(sizes[0]):
            for rest in cycle_profiles(sizes[1:]):
                yield (cluster,*rest)


def check_gcd_budget():
    summary={}
    for cost in range(8):
        histogram=Counter()
        extremal=set()
        for shape in shapes(cost):
            for cycles in cycle_profiles(shape):
                width=sum(max(c)//2 for c in cycles)
                assert width<=cost
                histogram[width]+=1
                if width==cost:
                    extremal.add(cycles)
        assert extremal=={((2,),)*cost}
        summary[cost]=dict(width_counts=sorted(histogram.items()))
    return summary


def determinant(matrix,p):
    a=[[x%p for x in row] for row in matrix]
    d=1
    for j,row in enumerate(a):
        first=next((k for k in range(j,len(a)) if a[k][j]),None)
        if first is None:return 0
        if first!=j:
            a[j],a[first]=a[first],a[j]
            d=-d
        d=d*a[j][j]%p
        inv=pow(a[j][j],-1,p)
        for target in a[j+1:]:
            scale=target[j]*inv
            for k in range(j,len(a)):
                target[k]=(target[k]-scale*a[j][k])%p
    return d%p


def falling(n,k):
    return prod(range(n-k+1,n+1))


def check_matrices(p):
    checks=0
    for width in range(1,7):
        for multiplicities in partitions(width):
            roots=tuple(p-1-3*j for j in range(len(multiplicities)))
            for first in range(-5,1-width+1):
                matrix=[]
                for q in range(first,first+width):
                    row=[]
                    for root,order in zip(roots,multiplicities):
                        for j in range(1,order+1):
                            coefficient=falling(q-1,j-1)*pow(factorial(j),-1,p)
                            row.append(coefficient*pow(root,q-j,p)%p)
                    matrix.append(row)
                assert determinant(matrix,p), (p,multiplicities,first)
                checks+=1
    for l in range(1,9):
        for j in range(l//2+1,p):
            assert 2*j>l
        assert 2*p-l*2>l
    return checks


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    assert not args.out.exists()
    t=sp.symbols("t")
    H=sp.Poly(t**6-4*t**5+6*t**4-7*t**3+4*t**2-2*t+1,t)
    m=1-t+t*t
    n=1-t*m
    b=1-(1-t)*n
    A=[1,t,0,n,b,0,0,0]
    B=[1,1-t,m,0,0,t*t*m,0,0]
    assert all((A[k]==0)!=(B[k]==0) or k==1 for k in range(1,6))
    assert all(sp.rem(sp.Poly(sum(A[i]*B[k-i] for i in range(k+1))-int(k<=6),t),H).is_zero
               for k in range(8))
    sums=[0]
    for k in range(1,8):
        sums.append(sp.expand(sum((-1)**(j-1)*A[j]*sums[k-j] for j in range(1,k))
                              +(-1)**(k-1)*k*A[k]))
    disc=int(sp.discriminant(H))
    assert disc==163*1279
    data={}
    for p,rho in ((163,-6),(1279,168)):
        h=sp.Poly(H.as_expr(),t,modulus=p)
        assert sp.gcd(h,h.diff())==sp.Poly(t-rho,t,modulus=p)
        assert h.eval(rho)==h.diff().eval(rho)==0 and h.diff().diff().eval(rho)!=0
        sk=[int(sp.diff(s,t).subs(t,rho))*pow(k,-1,p)%p for k,s in enumerate(sums) if k]
        assert sk[0]==1
        reciprocal=[1]
        for k in range(1,7):
            reciprocal.append(-sum(sk[i]*reciprocal[k-i] for i in range(1,k+1))%p)
        tail={e7:sp.Poly.from_list([*reciprocal,-e7],t,modulus=p) for e7 in range(1,p)}
        bad=[e7 for e7,K in tail.items() if sp.gcd(K,K.diff()).degree()!=0]
        data[p]=dict(collision_residue=rho%p,normalised_moments=sk,
                     hermite_sample_matrices_checked=check_matrices(p),
                     boundary_degree7_coefficients_before_constant=reciprocal,
                     boundary_nonzero_e7_with_repeated_root=bad,
                     boundary_squarefree_nonzero_e7_count=p-1-len(bad))
    result=dict(status="verified",sympy=sp.__version__,disc_H=disc,cycle_budget=check_gcd_budget(),
                primes=data,random_parameters=None,rounding=None)
    args.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(data,sort_keys=True))


if __name__=="__main__":
    main()
