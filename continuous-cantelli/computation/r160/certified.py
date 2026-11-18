#!/usr/bin/env python3
"""Directed-interval Bernoulli moment bracket for the infinite fat Cantor set."""
import json
from fractions import Fraction
from math import factorial
from pathlib import Path

import mpmath as mp

mp.iv.dps=60
iv=mp.iv
F=Fraction
N=100
M=42


def exact(x):
    return iv.mpf(x.numerator)/iv.mpf(x.denominator)


def bounds(a,b):
    return iv.mpf([exact(a).a,exact(b).b])


def make_moments():
    raw=[iv.mpf(1)]+[iv.mpf(0)]*(M+1)
    width=iv.mpf(1)
    keep=iv.mpf(1)
    for j in range(1,N+1):
        frac=iv.mpf(1)/(256*2**j)
        a=width*(1+frac)/4
        keep*=1-frac
        width*=((1-frac)/2)
        p=[iv.mpf(1)]
        for k in range(1,len(raw)):
            p.append(p[-1]*a*a/((2*k)*(2*k-1)))
        raw=[sum((raw[k-i]*p[i] for i in range(k+1)),iv.mpf(0))
             for k in range(len(raw))]
    return keep,[raw[k]*factorial(2*k) for k in range(len(raw))]


keep,mu=make_moments()
target=iv.sqrt(iv.mpf(8)/7)*iv.exp(-iv.mpf(2)/7)
mass_tail=iv.mpf(1)/(256*2**N)
width_tail=iv.mpf(1)/2**N
eseries=iv.mpf(1)/factorial(M+1)
# Uniform |p|≤2, |p''|≤16, |p''''|≤4000 for v≥1/8.
value_error=2*mass_tail+2*width_tail**2+2*eseries
curv_error=16*mass_tail+500*width_tail**2+48*eseries


def evaluate(v):
    coeff=iv.mpf(1)
    u,u2=iv.mpf(0),iv.mpf(0)
    for k in range(M+1):
        u+=coeff*mu[k]
        u2+=coeff*mu[k+1]
        coeff*=-iv.mpf(1)/(2*v*(k+1))
    factor=keep/iv.sqrt(2*iv.pi*v)
    value=u*factor+iv.mpf([-value_error.b,value_error.b])
    curv=(u2/v**2-u/v)*factor+iv.mpf([-curv_error.b,curv_error.b])
    return value,curv


lo,hi=F(1,8),F(1,6)
assert (evaluate(exact(lo))[0]-target).a>0
assert (evaluate(exact(hi))[0]-target).b<0
for _ in range(90):
    mid=(lo+hi)/2
    diff=evaluate(exact(mid))[0]-target
    if diff.a>0:
        lo=mid
    elif diff.b<0:
        hi=mid
    else:
        break
assert hi-lo<F(1,10**27)
assert (evaluate(exact(lo))[0]-target).a>0
assert (evaluate(exact(hi))[0]-target).b<0
curv=evaluate(bounds(lo,hi))[1]
assert curv.b<-2
ratio_value=2*target/(256*(-curv))
report=dict(N=N,taylor_order=M,mpmath=mp.__version__,precision=iv.dps,
            target=str(target),mass_interval=str(keep+iv.mpf([-mass_tail.b,0])),
            v_lo=str(lo),v_hi=str(hi),v_decimal_interval=str(bounds(lo,hi)),
            curvature_interval=str(curv),slope_interval=str(curv/2),
            root_difference_times_2_to_n_limit=str(ratio_value),
            critical_value_truncation_error=str(value_error))
Path("certificate.json").write_text(json.dumps(report,indent=2)+"\n")
for k,v in report.items():
    print(k,v)
