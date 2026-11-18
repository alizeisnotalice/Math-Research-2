"""Exact constants and Hermite algebra for the global Jensen sign examples."""
from fractions import Fraction as Q
from math import factorial

import sympy as S


y,z,r,g=S.symbols("y z r g",real=True)
P=1-y**2-g*y*(y**2-3)-r*(y**4-6*y**2+3)
assert S.expand(P+S.hermite_prob(2,y)+g*S.hermite_prob(3,y)+r*S.hermite_prob(4,y))==0
moment=lambda n: S.factorial2(n-1) if n%2==0 else 0
average=S.expand(P.subs(y,y+z))
assert S.expand(sum(average.coeff(y,k)*moment(k) for k in range(5))
                +z**2*(1+g*z+r*z**2))==0

# Beta=8, r=g=0. Explicit tilted second-moment coefficients.
q=Q(9,17)
c0=1-2*q+3*q*q
c2=(6*q-2)*Q(9,17**2)
c4=Q(81,17**4)
assert c0==Q(226,289)
small_A_half=Q(27,8)*c0          # actual limit is this divided by sqrt(17)
assert small_A_half**2<17
A=Q(1,64)
kappa=Q(4,153)-A/2
assert kappa>0
lower_e=Q(8,3)
uniform_upper=A*27/Q(16)*(c0+c2/(lower_e*kappa)+c4*(2/(lower_e*kappa))**2)
assert uniform_upper<Q(1,10)

# Omitted weighted tail for |G|>24, beta∈[1/64,512], |t|≤12.
# The universal |b|<2000 bound follows by optimizing each monomial term.
b_bound=1+65+Q(3)+Q(1000,6)+Q(130**2,12)
assert b_bound<2000
assert 513**3<12000**2
def upper_exp_negative(x):
    return 1/sum(x**j/factorial(j) for j in range(451))

min_linear_pref=2*Q(1,1024)*Q(15,16)/12000
ratio_tail=2000**2*upper_exp_negative(Q(144))/min_linear_pref
assert ratio_tail<Q(1,10**48)

phi2_remainder=Q(10,7)*Q(3,10)**7/factorial(9)
assert phi2_remainder<Q(9,10**10)
assert phi2_remainder/Q(7,20)<Q(3,10**9)
print("Hermite identities exact")
print("global J/(-s*m) upper bound for beta=8,A=1/64:",float(uniform_upper))
print("tail in normalized ratio over scan box:",float(ratio_tail))
print("small-z phi2 degree-six absolute error:",float(phi2_remainder))
