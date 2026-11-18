"""Reflection Laplace coefficient and two-moment bathtub checks."""
import json
import platform
import mpmath as mp
import sympy as sp

mp.mp.dps=48
T=mp.mpf('0.9')
B=mp.mpf('1.1')
A=mp.mpf('1.4')
p=lambda x:mp.exp(-x*x/(2*T))/mp.sqrt(2*mp.pi*T)
derivative=2*B*p(B)/T
def laplace_scaled(lam):
    return mp.quad(lambda z:mp.exp(-z)*lam*(p(B-z/lam)-p(B+z/lam)),[0,2,8,25,mp.inf])

# Use an independently feasible synthetic q to calibrate the two-moment cap.
q=lambda y:p(y)*mp.mpf('0.6')*(1-((2*y+A-B)/(A+B))**2)
integ=lambda fun,a,b:mp.quad(fun,[a,b])
M=integ(q,-A,B)
P=integ(lambda y:y*q(y),-A,B)
def residual(l,r):
    m=integ(p,-A,l)+integ(p,r,B)
    first=integ(lambda y:y*p(y),-A,l)+integ(lambda y:y*p(y),r,B)
    return m-M,first-P
l,r=mp.findroot(residual,(-mp.mpf('.7'),mp.mpf('.7')))
w=lambda y:(y+A)*(B-y)
bathtub=integ(lambda y:w(y)*p(y),-A,l)+integ(lambda y:w(y)*p(y),r,B)
actual=integ(lambda y:w(y)*q(y),-A,B)
assert bathtub<actual
u,t,b=sp.symbols('u t b',positive=True)
heat=lambda x:sp.exp(-x*x/(2*t))/sp.sqrt(2*sp.pi*t)
assert sp.simplify(sp.diff(heat(b-u)-heat(b+u),u).subs(u,0)-2*b*heat(b)/t)==0
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__,digits=mp.mp.dps,
                      synthetic_parameters=[str(A),str(B),str(T)],predicted_coefficient=mp.nstr(derivative,23),
                      scaled_laplace=[(lam,mp.nstr(laplace_scaled(lam),23)) for lam in (10,100,1000)],
                      cutpoints=[mp.nstr(l,20),mp.nstr(r,20)],bathtub=mp.nstr(bathtub,20),actual=mp.nstr(actual,20)),indent=2))
