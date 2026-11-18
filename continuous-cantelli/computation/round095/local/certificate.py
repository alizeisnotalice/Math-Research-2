"""Exact constants for the reflection crossover and integration cutoff."""
from fractions import Fraction as Q
from math import factorial

# For z in [0,1], the integrand defining M(1) is bounded below by e^-1.5.
neg=sum((-Q(3,2))**j/factorial(j) for j in range(6))
assert neg>Q(1,5)
example=Q(64**2,126)*(Q(1,5)-Q(1,127))
assert example>6

def e_upper_minus(x):
    return 1/sum(Q(x)**j/factorial(j) for j in range(151))

boundary_tail=max(49*e_upper_minus(48),e_upper_minus(98))
interior_tail=11*e_upper_minus(98)
small_sinhc=Q(1,4)**10/factorial(11)/(1-Q(1,16)/(12*13))
assert boundary_tail < Q(1,10**18)
assert interior_tail < Q(1,10**40)
assert small_sinhc < Q(1,10**13)
print("k^2 R at k=64,v=63 exceeds",float(example),">6")
print("endpoint normalized tail upper",float(boundary_tail))
print("interior normalized tail upper",float(interior_tail))
print("sinhc degree-eight real truncation upper",float(small_sinhc))

far_lower_loss=(Q(100,81)+Q(100,99)+1)/256
far_upper_excess=Q(1,99)
assert max(far_lower_loss,far_upper_excess)<Q(13,1000)
print("far weak wall |k^2 R-1| upper",float(max(far_lower_loss,far_upper_excess)))

import sympy as S

w,v2,h,tau,b,d=S.symbols("w v2 h tau b d", positive=True)
series=S.series(S.exp(-w*w/2)*(1+v2*w*w/6+v2*v2*w**4/120),w,0,6).removeO()
assert S.factorial(3)*series.coeff(w,2)==v2-3
assert S.factorial(5)*series.coeff(w,4)==v2*v2-10*v2+15
k_inverse_sq=h*(tau+h)/(tau*d*d)
scaled_v2=b*b*h/(tau*(tau+h))
physical=S.series(1+(scaled_v2-3)*k_inverse_sq+
    (scaled_v2**2-10*scaled_v2+15)*k_inverse_sq**2,h,0,3).removeO()
assert S.simplify(physical-(1-3*h/d**2+h**2*(b*b/(tau*tau*d*d)-3/(tau*d*d)+15/d**4)))==0
print("fixed-wall h and h^2 coefficients verified symbolically")
