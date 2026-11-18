"""Certified rational threshold boxes for one Gaussian-capped moment problem."""
from fractions import Fraction as Q
from pathlib import Path
import json
import math

import mpmath as mp

N=60


def plus(a,b): return a[0]+b[0],a[1]+b[1]
def minus(a,b): return a[0]-b[1],a[1]-b[0]
def times(a,b):
    c=[x*y for x in a for y in b]
    return min(c),max(c)


def eval_positive(x,integral):
    x=abs(x)
    terms=[((-1)**j)*x**(2*j)/(2**j*math.factorial(j)) for j in range(N+2)]
    if integral:
        terms=[x*t/(2*j+1) for j,t in enumerate(terms)]
    # Even and odd consecutive partial sums provide an exact bracket.
    return sum(terms[:N+2]),sum(terms[:N+1])


def ant(x):
    z=eval_positive(x,True)
    return z if x>=0 else (-z[1],-z[0])


def raw_density(box):
    near=min(abs(x) for x in box) if not box[0]<=0<=box[1] else Q(0)
    far=max(abs(x) for x in box)
    return eval_positive(far,False)[0],eval_positive(near,False)[1]


def approx_bbox(x,precision):
    n=int(mp.floor(x*precision))
    return Q(n-1,precision), Q(n+2,precision)


def atan(x,n=44):
    terms=[(-1)**j*x**(2*j+1)/(2*j+1) for j in range(n+2)]
    return sum(terms[:n+2]),sum(terms[:n+1])


def main():
    mp.mp.dps=70
    l,r=Q(-2),Q(3)
    mass,mu=Q(1,2),-Q(1,8)
    v=lambda q:mp.mpf(q.numerator)/q.denominator
    raw=lambda z:mp.exp(-z*z/2)
    A=lambda z:mp.sqrt(mp.pi/2)*mp.erf(z/mp.sqrt(2))
    def get(p):
        return (mp.findroot(lambda x:A(x)-A(-2)-p,(-1.5,-.1)),
                mp.findroot(lambda x:A(3)-A(x)-v(mass)+p,(.4,2.2)))
    def first(p):
        a,b=get(p)
        return raw(-2)-raw(a)+raw(b)-raw(3)
    p=mp.findroot(lambda x:first(x)-v(mu),(.2,.4))
    pbox=approx_bbox(p,10**14)
    assert 0<pbox[0]<pbox[1]<mass
    quantiles=[]
    for q in pbox:
        a,b=get(v(q))
        ax,bx=approx_bbox(a,10**20),approx_bbox(b,10**20)
        assert l<ax[0]<ax[1]<bx[0]<bx[1]<r
        assert minus(ant(ax[0]),ant(l))[1]<q<minus(ant(ax[1]),ant(l))[0]
        assert minus(ant(r),ant(bx[1]))[1]<mass-q<minus(ant(r),ant(bx[0]))[0]
        f=minus(plus(minus(raw_density((l,l)),raw_density(ax)),raw_density(bx)),raw_density((r,r)))
        quantiles.append((ax,bx,f))
    assert quantiles[0][2][0]>mu>quantiles[1][2][1]
    aa=(quantiles[0][0][0],quantiles[1][0][1])
    bb=(quantiles[0][1][0],quantiles[1][1][1])
    bound=(mass,mass)
    for sign,box in [(1,(l,l)),(-1,aa),(1,bb),(-1,(r,r))]:
        product=times(box,raw_density(box))
        bound=plus(bound,product) if sign>0 else minus(bound,product)
    pi=minus(times((Q(16),Q(16)),atan(Q(1,5))),times((Q(4),Q(4)),atan(Q(1,239))))
    scale=10**25
    sq=(Q(math.isqrt(math.floor(2*pi[0]*scale**2)),scale),
        Q(math.isqrt(math.floor(2*pi[1]*scale**2))+1,scale))
    norm=(bound[0]/sq[1],bound[1]/sq[0])
    assert norm[1]-norm[0]<Q(1,10**11)
    out=dict(lower_endpoint=str(l),upper_endpoint=str(r),
             raw_mass=str(mass),raw_mean_moment=str(mu),left_split=[str(z) for z in pbox],
             a=[str(z) for z in aa],b=[str(z) for z in bb],
             a_decimal=list(map(float,aa)),b_decimal=list(map(float,bb)),
             exact_normalized_second_interval=[str(z) for z in norm],
             normalized_second_decimal=list(map(float,norm)),
             second_interval_width=float(norm[1]-norm[0]),digits=mp.mp.dps)
    Path('exact_example.json').write_text(json.dumps(out,indent=2)+"\n")
    print(out['a_decimal'],out['b_decimal'],out['normalized_second_decimal'],out['second_interval_width'])

if __name__=='__main__': main()
