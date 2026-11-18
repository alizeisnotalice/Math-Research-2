"""Exact symbolic constants and direct quadrature for the cut-crossing bound."""
import json
import math
import sympy as sp
from scipy.integrate import quad
from scipy.special import ndtr
s=sp.Symbol('s', real=True)
assert sp.simplify(1+sp.integrate((1-s)/(1+s),(s,0,1))-2*sp.log(2))==0
density=lambda x: math.exp(-x*x/2)/math.sqrt(2*math.pi)
hermite=2*quad(lambda x: (1-x*x)*density(x),0,1)[0]+2*quad(lambda x: (x*x-1)*density(x),1,math.inf)[0]
assert abs(hermite-4*density(1))<1e-12
c=hermite*math.log(2)/2
checks=[]
for eps in (0.01,0.1,1.0):
    theta=math.atan(math.sqrt(eps))
    for z in (0,0.5,1.5):
        # Integrate in displacement from the threshold, in noise-scale units.
        scale=math.sin(theta)
        val=2*quad(lambda d: scale*density(z-scale*d)*ndtr(((math.cos(theta)-1)*z/scale)-math.cos(theta)*d),0,math.inf)[0]
        bound=theta/math.pi*math.exp(-z*z/2)
        assert val<=bound+1e-13
        if z==0: assert abs(val-bound)<1e-13
        checks.append([eps,z,val,bound])
print(json.dumps(dict(weighted_flip_lower=(1-c)/2,crossings=checks),indent=2))
