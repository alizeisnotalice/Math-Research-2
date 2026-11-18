"""Middle-third model, linear distance damping; integrate every equal-level gap exactly."""
import json, math, cmath
# Explicitly sum real parts around each parent midpoint, including the parent center phase.
# Parent centers have phase average (-1)^n times the product of prior cosines.
def coefficient(k):
    n=3**k; freq=2*math.pi*n; lam=float(n*n); c=-1.; a=b=0.
    for j in range(1,100):
        h=3.**(-j)/2
        z=complex(-lam,freq)*h
        H=(cmath.exp(z)-1)/complex(-lam,freq) if abs(z)>1e-5 else h*(1+z/2+z*z/6)
        angle=math.pi*((n%(2*3**j))/3**j)
        a+=2**j*c*(cmath.exp(-1j*angle)*H).real
        b+=2**j*(-math.expm1(-lam*h))/lam
        c*=math.cos(2*math.pi*((n%3**j)/3**j))
    return a/b
target=-math.prod(math.cos(2*math.pi/(3**j)) for j in range(1,40))
vals=[dict(k=k,ratio=coefficient(k)) for k in [2,4,6,8,10,12]]
assert abs(vals[-1]['ratio']-target)<.001
print(json.dumps(dict(limit=target,rows=vals),indent=2))
