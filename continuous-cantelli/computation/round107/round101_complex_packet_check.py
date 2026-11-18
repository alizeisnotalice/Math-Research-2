"""Compact-phase example and the translated complex diagonal, at 55 decimal digits."""
import json
import mpmath as m
import sympy as s
import platform

m.mp.dps=55
Z=lambda x:m.exp(-x*x/2)/m.sqrt(2*m.pi)
Phi=lambda x:(1+m.erf(x/m.sqrt(2)))/2
phase=lambda x:m.exp(-2*m.pi*1j*Phi(x))
mean=m.quad(lambda x:Z(x)*phase(x),[-m.inf,0,m.inf])
first=m.quad(lambda x:x*Z(x)*phase(x),[-m.inf,0,m.inf])
assert abs(mean)<m.mpf('1e-50') and first.imag>0
z,t,r,v,c=s.symbols('z t r v c',real=True, nonzero=True)
d=c+v/t
phase_full=(t+r)*(t+z)-(t+r)**2*d
phase_left=t*z-t*v+r*(z-2*v)-r*r*v/t
difference=s.simplify(phase_full-phase_left)
assert not difference.has(v,z)
print(json.dumps(dict(python=platform.python_version(),sympy=s.__version__,mpmath=m.__version__,digits=m.mp.dps,
                      phase_constant=str(difference),mean_modulus=m.nstr(abs(mean),12),
                      first_moment_imaginary=m.nstr(first.imag,30)),indent=2))
