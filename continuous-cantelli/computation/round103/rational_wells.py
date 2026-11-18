"""Exact rational certificate for two separated low-variance wells."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import json
import math

M=Q(9,4); e=Q(1,4); B=Q(1,64)
L=Q(1,16);R=Q(5,8)
S=10**12
n=256


def sqrt_bounds(x):
    t=math.isqrt((x*S*S).numerator//(x*S*S).denominator)
    return Q(t,S),Q(t+1,S)


def erf_exponential(x):
    assert 0<=x<=Q(45,2)
    term=Q(1);total=term
    for k in range(1,103):
        term*=-x/k
        total+=term
    return max(Q(0),total+term*(-x/103)),total


def atan(q):
    v=sum((-1)**j*q**(2*j+1)/(2*j+1) for j in range(20))
    return v,v+q**41/41


pi5,pi239=atan(Q(1,5)),atan(Q(1,239))
PI=(16*pi5[0]-4*pi239[1],16*pi5[1]-4*pi239[0])
normal=(sqrt_bounds(2*PI[0])[0],sqrt_bounds(2*PI[1])[1])


@lru_cache(maxsize=None)
def item(x,v):
    T=e+v
    a,b=sqrt_bounds((M+e)/T)
    l,u=erf_exponential((1+1/T)*x*x/2)
    return Q(math.floor(S*a*l),S),Q(math.ceil(S*b*u),S)


@lru_cache(maxsize=None)
def lower_envelope(x,U):
    p,q=item(x,Q(0)),item(x,U)
    return min(p[0],q[0]),min(p[1],q[1])


def half(use_well):
    low=high=Q(0)
    for i in range(3*n):
        a=Q(i,n);b=Q(i+1,n)
        U=B if use_well and L<=a<b<=R else M
        low+=lower_envelope(b,U)[0]/n
        high+=lower_envelope(a,U)[1]/n
    # The envelope is bounded by the narrow T=e branch; outside 3 use Mills.
    high+=item(Q(3),Q(0))[1]/(3*(1+1/e))
    return low,high


with_well=half(True)
bare=half(False)
joint_lower=2*with_well[0]/normal[1]
one_upper=(with_well[1]+bare[1])/normal[0]
mass_upper=Q(0)
for i in range(int(L*n),int(R*n)):
    left=Q(i,n)
    mass_upper+=erf_exponential(left*left/2)[1]/n
budget=2*(M-B)*mass_upper/normal[0]
assert joint_lower>Q(5,4)
assert one_upper<Q(97,100)
assert budget<Q(95,100)
out=dict(M=str(M),smoothing=str(e),y='0',wells=[['-5/8','-1/16'],['1/16','5/8']],well_upper=str(B),
         independent_cells_per_unit=n,ratio_two_wells_lower=str(joint_lower),
         ratio_one_well_upper=str(one_upper),mean_budget_upper=str(budget),
         decimal_ratios=[float(joint_lower),float(one_upper),float(budget)])
Path('rational_wells.json').write_text(json.dumps(out,indent=2)+"\n")
print(out['decimal_ratios'])
