"""Uniform exact integrated tail-well witness for all centers a >= 16.

Suppose f(a+y) <= 1 - (6/a^2)*(1-y^2)^2 for |y|<=1.
For y in [-1/4,1/4] and tau in [-1/8,1/8], t=a+tau,
the original positive kernel action is
  -(y-tau)^2/2 + (1-f(a+y))*t^2/2.
The rational lower bound is valid uniformly in a. Integrating the positive
x-core and J gives the same contribution as the Gaussian dual by Tonelli.
"""
from fractions import Fraction as F
import json
from math import factorial

minimum_center = F(16)
core_halfwidth = F(1,4)
interval_halfwidth = F(1,8)
depth = F(6)
squared_shape_lower = (1-core_halfwidth**2)**2
scaled_t_lower = 1-interval_halfwidth/minimum_center
action_lower = depth*scaled_t_lower**2*squared_shape_lower/2-(core_halfwidth+interval_halfwidth)**2/2
assert action_lower > F(5,2)
# sqrt(2*pi)<3 from pi<4.
density_normalization_lower = F(1,3)
exp_three_terms = 1+F(5,2)+F(5,2)**2/2
ratio_lower = 2*core_halfwidth*density_normalization_lower*exp_three_terms
assert ratio_lower == F(53,48)>1
# Refine the entire |y|<=1 core with one-sided endpoint bounds.
shallow_depth = F(8,5)
normalization_lower = F(39,100)
assert normalization_lower**2 * F(44,7) < 1
n = 64
scale = 2**24
units = 0
for i in range(1,n+1):
    u = F(i,n)
    e = shallow_depth*scaled_t_lower**2*(1-u**2)**2/2 - (u+interval_halfwidth)**2/2
    # Odd Taylor degree gives a global lower bound by Lagrange remainder.
    minimum_exp = sum(e**j/factorial(j) for j in range(10))
    units += (minimum_exp.numerator*scale)//minimum_exp.denominator
shallow_ratio_lower = normalization_lower * F(2,n) * F(units,scale)
assert shallow_ratio_lower == F(27467687481,26843545600) > 1
print(json.dumps({'centersAAtLeast':str(minimum_center),'depth':str(depth),'profileRadius':1,
   'intervalHalfWidth':str(interval_halfwidth),'rectangleHalfWidth':str(core_halfwidth),
   'actionLower':str(action_lower),'strictDualRatioLower':str(ratio_lower),
   'shallowDepth':str(shallow_depth),'shallowStrictDualRatioLower':str(shallow_ratio_lower)},indent=2))
