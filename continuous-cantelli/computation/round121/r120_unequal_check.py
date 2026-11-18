"""Exact algebra for the thresholds, and Gaussian tail-margin calibration."""
import json
import platform
import mpmath as mp
import sympy as sp

m,b,hp,hm,rhm,t,y2=sp.symbols("m b hp hm rhm t y2",nonzero=True)
cov=rhm-m*hp
lhs=hm+rhm/m-2*(1+t)
rhs=2*b*(m-b)+2*(b**2-y2)+cov/m
assert sp.simplify((lhs-rhs).subs(y2,1+t+b*m-(hp+hm)/2))==0
mp.mp.dps=50
rows=[]
for bb in (mp.mpf("0.1"),1,3):
    g=mp.exp(-bb**2/2)/mp.sqrt(2*mp.pi)
    tail=mp.erfc(bb/mp.sqrt(2))/2
    m=g/tail
    direct=mp.quad(lambda r:(r-bb)*mp.exp(-r**2/2)/mp.sqrt(2*mp.pi),[bb,mp.inf])
    assert abs(direct-(g-bb*tail))<mp.mpf("1e-47")
    mad=(mp.exp(-m**2/2)/mp.sqrt(2*mp.pi)-m*mp.erfc(m/mp.sqrt(2))/2)/tail
    rows.append({"forced_high_oscillation":mp.nstr(2*m*bb*(m-bb)/mad,18),
                 "radius":str(bb),"two_sided_budget":mp.nstr(2*bb*direct,18),
                 "weighted_conditional_margin":mp.nstr(2*bb*(m-bb),18)})
print(json.dumps(dict(python=platform.python_version(),sympy=sp.__version__,
                      mpmath=mp.__version__,digits=mp.mp.dps,
                      formula_exact=True,gaussian_margins=rows),indent=2))
