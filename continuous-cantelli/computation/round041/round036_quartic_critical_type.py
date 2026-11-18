"""Exact isolator and high-precision action/record checks for quartic threshold."""
from fractions import Fraction as F
from math import factorial
import json
import platform
import mpmath as mp

mp.mp.dps=65
def cubic(x):
    return x**3-4*x**2+6*x-2
lo,hi=F(1,3),F(1,2)
for _ in range(170):
    mid=(lo+hi)/2
    if cubic(mid)>0:
        hi=mid
    else:
        lo=mid
assert cubic(lo)<0<cubic(hi)
theta=mp.mpf((lo+hi).numerator)/(lo+hi).denominator/2
A=(1-theta)/mp.sqrt(theta*(mp.mpf('.5')-theta))
mu_over_l=4*mp.sqrt(mp.e/2)*(mp.mpf('.5')-theta)**mp.mpf('1.5')/(1-theta)
def qlog(x):
    return x*mp.log(x) if x else mp.mpf(0)
def phi(x):
    return x*mp.log(A*A)+x*mp.log(4*x)+qlog(mp.mpf('.5'))-qlog(mp.mpf('.5')-x)+qlog(2-2*x)-2*mp.log(2)
assert abs(phi(theta))<mp.mpf('1e-47')
# Roots / thresholds are checked with true rational arithmetic by squaring
#  A=(1-theta)/sqrt(theta(1/2-theta)) in this interval; denominator
# is decreasing, and numerator is decreasing.
lb,ub=F(385,100),F(386,100)
assert (1-hi)**2>lb**2*lo*(F(1,2)-lo)
assert (1-lo)**2<ub**2*hi*(F(1,2)-hi)
# Check the Borel multiplication and its inverse without approximating series.
src=[F((-1)**j*(j+2),3**j) for j in range(25)]
c=F(7,5)
def gaussian_convolve(coefficients,sign):
    return [sum((sign*c/2)**i/factorial(i)*coefficients[k-2*i] for i in range(k//2+1)) for k in range(len(coefficients))]
assert gaussian_convolve(gaussian_convolve(src,1),-1)==src
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,
  exact_theta_bracket=[str(lo),str(hi)],theta=mp.nstr(theta,32),
  critical_type_over_sqrt_c=mp.nstr(A,32),verified_type_bracket=["3.85","3.86"],
  mu_over_lambda=mp.nstr(mu_over_l,28),
  correction_over_lambda=mp.nstr((1-mu_over_l)/mp.sqrt(2*theta),28),
  tangent_action=mp.nstr(phi(theta),10),exact_inverse_terms=len(src)),indent=2))
