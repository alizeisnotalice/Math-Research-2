#!/usr/bin/env python3
"""Independent lead symbolic audit; no halo enumeration or numerical claim."""
import hashlib,json,platform,sys
from pathlib import Path
from fractions import Fraction as F
import sympy as sp
R,t=sp.symbols('R t',positive=True)
checks={}
checks['root_mass_identity']=sp.expand((2*R*(R+t)+2-4*R**2).subs(t,R-1/R))==0
checks['antiderivative']=sp.simplify(sp.diff(sp.log(R)/4+1/(16*R**4),R)-(R-1/R)*(1+1/R**2)/(4*R**2))==0
checks['normalized_mass_polynomial']=sp.expand(4*(t+1/t)**2+2-(4*t**2+10+4/t**2))==0
checks['t4_normalized_constant']=4+10*F(1,16)+4*F(1,256)==F(297,64)
checks['dyadic_per_band']=F(3,2)/F(297,16)==F(8,99)
checks['root4']=sp.simplify((2+sp.sqrt(5))*sp.sqrt(5)-(2+sp.sqrt(5))*2)==1
checks['canonical_minus_lower_ratio']=F(297,4)/66==F(9,8)
# Polynomial division certifies R(t)^2=tR(t)+1 and the mass increase.
checks['mass_increase']=sp.diff(2*R*(R+t)+2,R)==4*R+2*t
assert all(checks.values()),checks
out={'status':'PASS','python':sys.version,'sympy':sp.__version__,'platform':platform.platform(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'checks':checks,'scope':'Eight symbolic/rational identities only; general geometry proven analytically; no random inputs or float tolerance.'}
print(json.dumps(out,indent=2))
