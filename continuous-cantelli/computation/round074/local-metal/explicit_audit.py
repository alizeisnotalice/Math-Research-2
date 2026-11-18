"""Check the compressed level formula against all finite gaps separately."""
import cmath
import json
import math
from pathlib import Path

from scipy.integrate import quad

from verify import level_integral


def verify(lam,power,d,n,depth):
    intervals=[(0.0,1.0)]
    observed=0j
    a=2*math.pi**2*n*n*d
    for _ in range(depth):
        new=[]
        for l,r in intervals:
            h=(r-l)*lam
            x0,x1=l+h,r-h
            width=x1-x0
            half=width/2
            for endpoint,sgn in ((x0,1),(x1,-1)):
                def term(t):
                    u=half*t*t
                    return cmath.exp(2j*math.pi*n*(endpoint+sgn*u))*math.expm1(-a*u**power)*2*half*t
                observed+=complex(quad(lambda t:term(t).real,0,1,epsabs=1e-14)[0],
                                  quad(lambda t:term(t).imag,0,1,epsabs=1e-14)[0])
            new.extend(((l,x0),(x1,r)))
        intervals=new

    grouped=0.0
    phase_prod=1.0
    for j in range(depth):
        length=(1-2*lam)*lam**j
        omega=math.pi*n*length
        s=a*(length/2)**power
        grouped+=(-1)**n*(2**j)*length*phase_prod*level_integral(power,omega,omega,s)
        phase_prod*=math.cos(math.pi*n*(1-lam)*lam**j)
    discrepancy=abs(observed-grouped)
    assert discrepancy<1e-12
    return {"lambda":lam,"p":power,"d":d,"n":n,"level":depth,
            "gaps":2**depth-1,"explicit":[observed.real,observed.imag],
            "grouped":grouped,"absolute_difference":discrepancy}


cases=[verify(1/3,.75,2,27,8),verify(.41,1.5,.125,43,8)]
Path("explicit-audit.json").write_text(json.dumps(cases,indent=2)+"\n")
print(cases)
