#!/usr/bin/env python3
"""Exact rank audit of the finite moment windows, with jet-derived moments."""

import argparse
import json
from pathlib import Path


def mul(a,b,p):
    return a[0]*b[0]%p,(a[0]*b[1]+a[1]*b[0])%p


def add(a,b,p):
    return tuple((x+y)%p for x,y in zip(a,b))


def scale(m,a,p):
    return tuple(m*x%p for x in a)


def moments(p,t0):
    t=t0%p,1
    one=(1,0)
    zero=(0,0)
    t2=mul(t,t,p)
    t3=mul(t2,t,p)
    n=add(add(add(one,scale(-1,t,p),p),t2,p),scale(-1,t3,p),p)
    b=add(one,scale(-1,mul(add(one,scale(-1,t,p),p),n,p),p),p)
    a=(one,t,zero,n,b,zero,zero,zero)
    ps=[zero]
    for k in range(1,8):
        out=scale((-1)**(k-1)*k,a[k],p)
        for i in range(1,k):
            out=add(out,scale((-1)**(i-1),mul(a[i],ps[k-i],p),p),p)
        ps.append(out)
    return [ps[k][1]*pow(k,-1,p)%p for k in range(1,8)]


def rank(rows,p):
    a=[row[:] for row in rows]
    s=0
    for k in range(len(a[0])):
        slot=next((i for i in range(s,len(a)) if a[i][k]),None)
        if slot is None:continue
        a[s],a[slot]=a[slot],a[s]
        inverse=pow(a[s][k],-1,p)
        for i in range(s+1,len(a)):
            fac=a[i][k]*inverse%p
            for j in range(k,len(a[i])):
                a[i][j]=(a[i][j]-fac*a[s][j])%p
        s+=1
    return s


def audit(p,rho):
    u=moments(p,rho)
    assert u[0]==1
    coeffs=[1]
    for k in range(1,7):
        coeffs.append(-sum(coeffs[k-i]*u[i] for i in range(1,k+1))%p)
    K=[0]+coeffs[::-1]
    assert K[0]==0 and K[-1]==1
    assert coeffs=={163:[1,6,121,106,68,115,117],1279:[1,1111,1202,215,1193,762,576]}[p]
    results=[]
    for s in range(6,13):
        for m in range(13):
            matrix=[[u[i+k] if 0<=i+k<7 else 0 for i in range(m+1)]
                    for k in range(-s,7-m)]
            if m<=s:
                assert [0]*m+[1] in matrix
                continue
            degree_q=m-s-1
            degree_r=m-7
            expected_free=degree_q+degree_r+1
            assert rank([r[:-1] for r in matrix],p)==m-expected_free
            apart=(s-6)*[0]+K
            lift=[0]*degree_q+apart
            assert len(lift)==m+1
            def vanishes(coeff):
                return all(sum(a*b for a,b in zip(coeff,row))%p==0 for row in matrix)
            assert vanishes(lift)
            assert all(vanishes([0]*j+apart) for j in range(degree_q))
            assert all(vanishes([0]*j+[1]) for j in range(degree_r+1))
            results.append(dict(m=m,s=s,free_coefficients=expected_free,
                                nonzero_constant_family_size=p**(expected_free-1)*(p-1)))
    return dict(moments=u,K_ascending=K,per_window=results)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out', required=True, type=Path)
    out=parser.parse_args().out
    assert not out.exists()
    payload=dict(status='verified',instances={p:audit(p,rho) for p,rho in ((163,-6),(1279,168))})
    out.write_text(json.dumps(payload,sort_keys=True,indent=2)+'\n')
    print(json.dumps({p:len(e['per_window']) for p,e in payload['instances'].items()}))


if __name__=='__main__':
    main()
