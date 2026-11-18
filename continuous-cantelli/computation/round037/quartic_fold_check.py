"""Exact fold jets and numerical matched sqrt(N) amplification in the rational model."""
import json
import platform
import mpmath as mp
import sympy as sp

u,e,a,h2,b,w,y=sp.symbols('u e a h2 b w y')
F=u**2*(1+a*u+h2*u**2)+e*u
hc=-e/2-3*a*e**2/8
assert sp.series(sp.diff(F,u).subs(u,hc),e,0,3).removeO().expand()==0
Fc=sp.series(F.subs(u,hc),e,0,4).removeO().expand()
assert sp.expand(Fc+e**2/4+a*e**3/8)==0
assert sp.series(-e/Fc,e,0,1).removeO().expand()==4/e-2*a
quadratic=(b-a)*u**2+(1+a*w)*u-w
assert sp.expand(sp.discriminant(quadratic,u)-(1+(4*b-2*a)*w+a*a*w*w))==0
mp.mp.dps=55
d=mp.mpf('0.8')
aa=mp.mpf('1.3')
A=aa/mp.sqrt(d)
def compute(N):
    mode=mp.sqrt(d*N)
    def logweight(offset):
        r=mode+offset
        return N*mp.log(r/mode)-(r*r-mode*mode)/(2*d)
    def ref(offset):
        s=mode+offset
        # b=-s; compute the larger absolute discriminant factor.
        kappa=2*s+aa+2*mp.sqrt(s*(s+aa))
        return mp.exp(logweight(offset)+N*mp.log(kappa/(4*s))-A*mp.sqrt(N)/2)
    den=mp.quad(lambda z:0 if z<=-mode else mp.exp(logweight(z)),[-mode,-1,0,1,10,mp.inf])
    ref_mass=mp.quad(ref,[0,1])/den
    t0=mp.sqrt(4*N-16)
    def second(z,tilt):
        if z<=-t0:
            return 0
        return mp.exp((2*N-8)*mp.log1p(z/t0)-(2*t0*z+z*z)/4 +(A*z/4 if tilt else 0))
    panels=[-t0,-1,0,1,10,mp.inf]
    tilt=mp.exp(A*t0/4-A*mp.sqrt(N)/2)*mp.quad(lambda z:second(z,True),panels)/mp.quad(lambda z:second(z,False),panels)
    return dict(N=N,normalized_rational_fold_unit_window=mp.nstr(ref_mass,18),normalized_Hermite_bound_tilt=mp.nstr(tilt,18))
print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'digits':mp.mp.dps,'a':str(aa),'d':str(d),'critical_value':str(Fc),
    'exact_identities':4,'rows':[compute(n) for n in (32,128,512,2048)]},indent=2))
