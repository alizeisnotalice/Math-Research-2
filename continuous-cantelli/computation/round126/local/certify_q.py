#!/usr/bin/env python3
"""Interval Riemann enclosure for one Gaussian-smoothed Q witness."""
import json
import struct
from fractions import Fraction
from pathlib import Path

from mpmath import iv

iv.dps=45
norm=1/iv.sqrt(2*iv.pi)


def f32q(x):
    return Fraction.from_float(struct.unpack("=f",struct.pack("=f",x))[0])


def exact(x):
    return iv.mpf(x.numerator)/iv.mpf(x.denominator)


def hull(left,right):
    return iv.mpf([exact(left).a,exact(right).b])


def density(x,v):
    return norm/iv.sqrt(v)*iv.exp(-(x**2)/(2*v))


def integral(l,r,n,fn):
    if l>=r:
        return iv.mpf(0)
    h=(r-l)/n
    total=iv.mpf(0)
    for j in range(n):
        x=hull(l+j*h,l+(j+1)*h)
        total+=fn(x)*exact(h)
    return total


data=json.loads(Path("scan-main.json").read_text())
profiles=json.loads(Path("candidates.json").read_text())["candidates"]
case=len(profiles)-1
profile=profiles[case]
selected=data["byProfile"][case][6]["best"]
n=64
t,y,eps=map(f32q,(selected["p"],selected["q"],selected["aux"]))
assert t>=1 and eps>0
domain=f32q(profile["R"])
dx=f32q(profile["dx"])
b=list(map(f32q,profile["offset"]))
total=iv.mpf(0)
for k,(b0,b1) in enumerate(zip(b,b[1:])):
    if b0==b1==0:
        continue
    l=-domain+k*dx
    r=l+dx
    y_iv=exact(y)
    total+=integral(l,r,n,lambda x:density(x,1)*density(y_iv-x,exact(t-1+eps)))
    slope=(b1-b0)/dx
    if slope==0:
        if b0>t-1:
            continue
    else:
        cutoff=l+(t-1-b0)/slope
        if slope>0:
            r=min(r,cutoff)
        else:
            l=max(l,cutoff)
    tp=exact(t-1+eps-b0)
    left=exact(-domain+k*dx)
    sr=exact(slope)
    total-=integral(l,r,n,lambda x:density(x,1)*density(y_iv-x,tp-sr*(x-left)))
assert total.b<iv.mpf("-.07")
print("profile",case,"t",float(t),"y",float(y),"epsilon",float(eps))
print("interval smoothed Q",total)
print("strict upper bound < -0.07")

for mode in (4,5):
    tau=f32q(data["byProfile"][case][mode]["best"]["p"])
    coefficient=exact(-tau*tau/2)
    tau_iv=exact(tau)
    residual=iv.mpf(0)
    for k,(b0,b1) in enumerate(zip(b,b[1:])):
        if b0==b1==0:
            continue
        l=-domain+k*dx
        r=l+dx
        offset, left, slope=exact(b0),exact(l),exact((b1-b0)/dx)
        residual+=integral(l,r,n,lambda x:density(x-tau_iv,1)
                           *(iv.exp(coefficient*(offset+slope*(x-left)))-1))
    assert residual.b<iv.mpf("-.01")
    print("held-out real",float(tau),residual)
    print("strict upper bound < -0.01")
