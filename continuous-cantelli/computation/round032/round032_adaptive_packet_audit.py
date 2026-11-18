"""Crosscheck the compact-spectrum packet at both time signs by independent quadratures."""
import cmath
import json
import math
import platform
import scipy
from scipy.integrate import quad
from scipy.stats import irwinhall

def psi(z, n):
    y=z/(4*n)
    return (cmath.sin(y)/y)**(2*n) if y else 1.

def row(n,z,d):
    B=irwinhall(2*n)
    def spec(k):
        return (2*n)*B.pdf(2*n*k+n)*math.exp(d*(1-k*k)/2)*math.cos(k*(z+d))
    # B-spline support scaled from [0,2n] to [-1/2,1/2].
    freq = quad(spec,-.5,.5,points=[(i-n)/(2*n) for i in range(1,2*n)],epsabs=1e-12)[0]
    def cond(g):
        return (math.exp(d/2-g*g/2)/math.sqrt(2*math.pi)
                * psi(z+d+cmath.sqrt(d)*g,n).real)
    phys = quad(cond,-15,15,epsabs=1e-12,epsrel=1e-12)[0]
    assert abs(phys-freq)<1e-9
    return dict(N=n,z=z,D=d,frequency=freq,physical=phys,error=abs(freq-phys))

print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,
    rows=[row(n,z,d) for n,z,d in ((1,0,0),(1,3,-5),(2,-6,2),(4,13,-10),
                                    (7,-22,7),(12,60,-23))]),indent=2))
