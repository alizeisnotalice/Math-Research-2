#!/usr/bin/env python3
"""Exact all-parameter identities and rational constants for r>=72."""
import hashlib
import json
from pathlib import Path
import platform
import sympy as S

s,y,r = S.symbols('s y r')
D = 1-2*s+3*s*s
T = (1+s**3*y)/((1+s)*D)
F = (1-s*s*y*y)/(1+s)-T*(1-(r+3)*s*s*y-(r+2)*s**3*y)/(1+s)**2+T*T*(1-3*s)
P = 9*s**6+6*s**5-5*s**4+4*s**3-3*s*s-2*s+1
A0 = 3*s**3+s*s-s+1
C0 = 3-3*s+3*s*s-3*s**3
A1 = 3*s**6+s**5-s**4+s**3
C1 = -6*s**6-3*s**5+2*s**4-s**3-3*s*s+2*s-1
Delta = s*s*y*(r*A0+C0+y*(r*A1+C1))/((1+s)**3*D*D)
assert S.cancel(F-P/((1+s)**3*D*D)-Delta) == 0
W4 = 1-2*s+3*s*s-4*s**3+5*s**4
L = 2*s**4-2*s**3+2*s*s-2*s+1
b4 = (1+s**5*y)/(1+s)-T*W4
assert S.cancel(b4-(-s**3*(5*s-4)-s**3*y*L)/((1+s)*D)) == 0
assert S.expand(L-((1-s)**2+s*s*(1-2*s+2*s*s))) == 0
assert S.expand(L-(1-2*s*(1-s)*(1+s*s))) == 0
lo,hi = S.Rational(13,20),S.Rational(17,20)
constants = {
 'lower_s_margin': S.Rational(1,20)*(1-lo**10)/(1+lo)-lo**9,
 'upper_s_margin': S.Rational(1,4)-hi**9,
 'P_prime_at_lo': S.diff(P,s).subs(s,lo),
 'P_at_lo':P.subs(s,lo),
 'P_second_positive_part_1':(270*s**4-60*s*s).subs(s,lo),
 'P_second_positive_part_2':(120*s**3+24*s-6).subs(s,lo),
 'denominator_lower':((1+s)**3*D*D).subs(s,lo),
 'F_infty_lower':P.subs(s,lo)/((1+s)**3*D*D).subs(s,hi),
 'tail_at_72':148*hi**74,
 'tail_ratio_at_4':S.Rational(14,12)*hi,
 'A0_absolute_bound':3*hi**3+hi**2+hi+1,
 'A1_absolute_bound':3*hi**6+hi**5+hi**4+hi**3,
 'C1_absolute_bound':6*hi**6+3*hi**5+2*hi**4+hi**3+3*hi**2+2*hi+1,
}
assert all(constants[k]>0 for k in ['lower_s_margin','upper_s_margin','P_prime_at_lo',
                                  'P_at_lo','P_second_positive_part_1','P_second_positive_part_2'])
assert constants['denominator_lower']>4
assert constants['F_infty_lower']==S.Rational(17168820,17453453557)
assert constants['tail_at_72']<constants['F_infty_lower']
assert constants['tail_ratio_at_4']<1
assert constants['A0_absolute_bound']<5
assert constants['A1_absolute_bound']<3
assert constants['C1_absolute_bound']<11
print(json.dumps({'python':platform.python_version(),'sympy':S.__version__,
                  'full_tail_identity':True,'b_r_plus_4_identity':True,
                  'positive_L_decomposition':True,
                  'constants':{k:str(v) for k,v in constants.items()},
                  'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
