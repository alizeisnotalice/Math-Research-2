"""Deterministic calibration with v=1-cos(x); this profile is nonstationary."""
import json
import platform
import scipy
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import ive
from math import cos, pi

def moment(t, j):
    z=t*t/2
    return ive(abs(j),z)/ive(0,z)

def mean_v(t):
    return 1-moment(t,1)

def energy_weight(t):
    return t*((1+moment(t,2))/2-moment(t,1)**2)

rows=[]
for a,b in ((1.,4.),(4.,16.),(16.,32.)):
    reg,err=quad(energy_weight,a,b,epsabs=1e-12,epsrel=1e-12)
    exact=mean_v(a)-mean_v(b)
    rows.append(dict(a=a,b=b,quadrature=reg,endpoint_difference=exact,discrepancy=abs(reg-exact),reported_error=err))
phase=minimize_scalar(lambda x:x*x+pi*cos(x),bounds=(0,pi),method="bounded",options={"xatol":1e-13})
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,profile="1-cos(x), calibration only",dissipation=rows,phase_minimizer=float(phase.x),phase_minimum=float(phase.fun),exact_minimum=pi**2/4),indent=2))
