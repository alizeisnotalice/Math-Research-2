"""Deterministic audit; floating values are not global proof certificates."""
import json, sys, platform
from fractions import Fraction
from math import factorial
import numpy as np
import sympy as sp
import scipy
from scipy.integrate import quad
from pathlib import Path
k=3.0
r=np.exp(-k*k/2)
B=4*k*k*r/(1-r*r)**2
m=1.0
phi=lambda x: np.exp(-x*x/2)/np.sqrt(2*np.pi)
a=lambda x: m+B*(np.cos(k*x)-r)
def E(f):
    v,e=quad(lambda x: f(x)*phi(x),-12,12,epsabs=2e-12,epsrel=2e-12,limit=400)
    return v,e
vals={}
for name,fun in [('mean',a),('fourth',lambda x:x**4+6*x*x*a(x)+3*a(x)**2),('sixth',lambda x:x**6+15*x**4*a(x)+45*x*x*a(x)**2+15*a(x)**3)]:
    vals[name]=E(fun)
transforms=[]
for t in [0.5,1.0,2.0,4.0]:
    v,e=E(lambda z:np.exp(t*t*(a(z+t)-m)/2))
    transforms.append(dict(t=t,shifted_mgf=v,quad_error=e))
M=sp.Matrix([[sp.cos(sp.pi*j*l/4)-1 for l in [1,2,3]] for j in [1,2,4]])
det=sp.simplify(M.det())
assert det!=0
# Exact Gaussian conditional-moment coefficients by binomial expansion.
x,A,t=sp.symbols('x A t')
polys={}
for n in [2,4,6]:
    expr=sum(sp.binomial(n,2*j)*x**(n-2*j)*A**j*sp.factorial2(2*j-1) for j in range(n//2+1))
    polys[str(n)]=str(sp.expand(expr))
assert abs(vals['mean'][0]-m)<1e-10
assert abs(vals['fourth'][0]-3*(1+m)**2)<1e-9
assert abs(vals['sixth'][0]-15*(1+m)**3)>1e-2
rr,kk=sp.symbols('r k', positive=True)
BB=4*kk**2*rr/(1-rr**2)**2
moment_identity=sp.factor(BB**2*(1-rr**2)**2/2-2*BB*kk**2*rr)
assert moment_identity==0
exp_partial=sum(Fraction(9,2)**j/factorial(j) for j in range(8))
assert exp_partial>80
out=dict(exact_fourth_identity=str(moment_identity),exp_4_5_lower=str(exp_partial),runtime=dict(python=platform.python_version(),executable=sys.executable,numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),seed=None,method='deterministic quad on [-12,12]; no interval certification',parameters=dict(k=k,B=B,m=m),range_a=[1-B*(1+r),1+B*(1-r)],moments=vals,fourth_residual=vals['fourth'][0]-12,sixth_residual=vals['sixth'][0]-120,transforms=transforms,cosine_matrix=str(M),exact_cosine_determinant=str(det),conditional_moment_polynomials=polys)
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
