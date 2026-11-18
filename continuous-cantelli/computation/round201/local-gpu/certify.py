#!/usr/bin/env python3
"""Interval mixed-sign check and rational uniform negative-tail criterion."""
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp

F=Fraction
iv=mp.iv
iv.dps=40
steps=768


def e(x):
    return iv.mpf(x.numerator)/iv.mpf(x.denominator)


def box(a,b):
    return iv.mpf([e(a).a,e(b).b])


def J(z,k):
    a=k+.5
    k=iv.mpf(k)
    c=1-k
    aa=iv.mpf(a)
    pref1=4*iv.mpf("0.5")**aa
    pref2=4*iv.mpf("0.5")**c
    z2=z*z
    acc=iv.mpf(0)
    for i in range(steps):
        t=box(F(i,steps),F(i+1,steps))
        q1=t**4/2
        q2=1-q1
        f1=pref1*t**int(4*a-1)*q2**(-k)
        f2=pref2*t**int(4*float(c.a)-1)*q2**(aa-1)
        acc+=f1*(1-z2*q1)*iv.exp(-z2*q1/2)
        acc+=f2*(1-z2*q2)*iv.exp(-z2*q2/2)
    return acc/(2*k*steps)


def certificate(m,D,k):
    rho=F(1,m*m)
    a=[F(9,10)+i*(F(1,10)-rho)/(m-2) for i in range(m-1)]
    z=J(box(F(0),D*rho),k)
    # Enclose full conditional Cantor support in each group of cylinders.
    bins=3 if m==64 else 4
    for j in range(0,m-1,bins):
        h=min(j+bins,m-1)
        z+=(h-j)*J(box(D*a[j],D*(a[h-1]+rho)),k)
    positive_branch_sum=z/m
    power=1+2*k
    r=iv.mpf(m)**(-int(4*k+1))
    tail=iv.mpf(100)*e(F(9*D,10))**(-power)*r/(1-r)
    return positive_branch_sum+iv.mpf([-tail.b,0]),positive_branch_sum,tail


# e>8/3 gives exp(-z^2/2) < eps for z>=4.5. The same
# monotone-weight comparison proves J<0 on all outer shells.
eps=F(3,8)**10
k_half_right=2*(1-eps)/225-4*eps
assert 35**4*2**3>225**3
k_quarter_right=(F(4,3)*F(4,11)+F(4,5)*F(16,121))/35 - 4*eps
assert k_half_right*255>2
assert k_quarter_right*255>4

rows=[]
for m,D in ((64,10),(128,15)):
    for k in (.25,.5):
        sk,inner,tail=certificate(m,D,k)
        if k==.25:
            assert sk.b<0
        else:
            assert sk.a>0
        rows.append(dict(m=m,D=D,k=k,all_scales=str(sk),inner=str(inner),outer_abs_bound=str(tail)))
out=dict(steps_per_half=steps,interval_digits=iv.dps,mpmath=mp.__version__,
         minimum_branch_count_uniform_negative=256,
         right_tail_upper_magnitudes_rational={"k1_4":str(k_quarter_right),"k1_2":str(k_half_right)},
         mixed_sign_intervals=rows)
Path("certificate.json").write_text(json.dumps(out,indent=2)+"\n")
for row in rows:
    print(row,flush=True)
