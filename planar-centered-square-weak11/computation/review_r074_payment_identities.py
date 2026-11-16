"""Symbolic branch identities for the independently read continuum payment."""
from pathlib import Path
import hashlib, json, os, platform, sys, time
import sympy as s
start=time.monotonic()
u,v,A,B,m,M,e,k=s.symbols('u v A B m M e k', real=True)
d=u-v
identities={
 'both_overlap_lengths_above_d':u*u-v*v-(A*B-(A-d)*(B-d))-d*(2*u-A-B),
 'one_overlap_length_at_most_d':u*u-v*v-d*u-d*v,
 'same_outside_mass_annulus_payment':(m+M)/4-M/4-m/4,
 'container_minus_corners':(2-e)**2-4*(s.Rational(1,2)-e)**2-(3-3*e**2),
 'even_grid_corner_budget':2*k*k/(2*k)**2-s.Rational(1,2),
 'odd_grid_corner_deficit':s.Rational(1,2)-2*k*k/(2*k+1)**2-(4*k+1)/(2*(2*k+1)**2),
}
zeros={name:str(s.factor(expr)) for name,expr in identities.items()}
assert set(zeros.values())=={'0'}
result={'status':'PASS','symbolic_zero_identities':zeros,
 'inequality_inputs':'0<=v<=u,0<=A,B<=u. If A<=d or B<=d, AB<=du. Otherwise annulus gap=d(2u-A-B)>=0. Geometric overlap shrink and localization are analytic proof review, not established by these identities.',
 'scope':'six algebraic identities; author 255/4335 regressions not rerun',
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'runtime':{'python':platform.python_version(),'sympy':s.__version__,'pid':os.getpid(),'seconds':time.monotonic()-start}}
Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
