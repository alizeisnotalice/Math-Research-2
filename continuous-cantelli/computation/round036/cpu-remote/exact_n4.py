"""Exact rational Sturm checks for an explicit n=4 Gaussian packet filter."""
import json
from pathlib import Path
import platform

import sympy as s

x = s.Symbol("x", real=True)
v = s.Symbol("v", real=True)
g = s.sqrt(3)/2
I = s.I
A = s.Rational(3, 8)
wstar = -s.Rational(1, 4)+I*g/2
delta = s.Rational(79, 100)
target_D = s.Rational(149, 1000)
cutoff=s.Rational(199,200)
outer=24*g


def rational_real(z):
    denominator = 1+delta*z**4/2
    numerator = (I*z-z**4/4)*denominator + delta*(I*z-z**4/2-wstar)**2/2
    return s.expand(s.re(s.expand(numerator*s.conjugate(denominator)))), s.expand(
        s.re(s.expand(denominator*s.conjugate(denominator))))


def strict(poly, variable, a, b, sign):
    p = s.Poly(poly, variable, domain=s.QQ)
    assert p.eval(0)*sign>0
    roots=p.count_roots(a,b)
    assert roots==0
    return dict(degree=p.degree(),sign_at_zero=int(s.sign(p.eval(0))),roots=int(roots))


real = g*x
den = 1+delta*real**4/2
real_num = -real**4/4 + delta*(-real**4/8+s.Rational(1,32)-(real-g/2)**2/2)
# D/z² exceeds A+Re G for every real z≠0.
real_poly = s.expand(target_D*den-real**2*(A*den+real_num))
real_check=strict(real_poly,x,-s.oo,s.oo,1)
num,den = rational_real(g*x+I/2)
horz = s.Poly(s.expand(num+A*den),x,domain=s.QQ)
assert horz.eval(1)==0 and horz.diff().eval(1)==0
quotient,remainder=s.div(horz,s.Poly((x-1)**2,x,domain=s.QQ))
assert remainder.is_zero
horizontal_check=strict(quotient.as_expr(),x,-s.oo,s.oo,-1)
vertical_checks=[]
for sign in (-1,1):
    num,den = rational_real(sign*outer+I*v)
    vertical_checks.append(strict(s.expand(num+A*den),v,0,s.Rational(1,2),-1))
# In 0≤Im z≤1/2, Re(z^4) ≥ -8(Im z)^4 ≥ -1/2.
assert 1-delta/4>0
# q=1/4, truncation outside |tau|<=cutoff lambda.
assert cutoff>delta/4
frequency_tail_action = cutoff/4-cutoff**2/(2*delta)
assert frequency_tail_action < -A
# Use the R027 large-input bound with a rectangular source cutoff L=2sqrt(3).
outer_source_action = -(1-cutoff)**2*outer**4/(2*4)+delta/32
assert outer_source_action < -A
# Implied original coefficient improvement normalized to the round-27 exact D*.
assert target_D>0 and target_D**2 < s.Rational(1,32)
out = dict(delta=str(delta),frequency_cutoff=str(cutoff),real_source_cutoff=str(outer),explicit_D=str(target_D),old_D="1/(4 sqrt(2))",
           relative_coefficient=str(s.simplify(target_D*4*s.sqrt(2))),
           real_action_gap_check=real_check,horizontal_quotient_check=horizontal_check,
           vertical_checks=vertical_checks,frequency_tail_action=str(frequency_tail_action),
           outer_source_action=str(outer_source_action),python=platform.python_version(),sympy=s.__version__)
Path('exact_n4.json').write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
