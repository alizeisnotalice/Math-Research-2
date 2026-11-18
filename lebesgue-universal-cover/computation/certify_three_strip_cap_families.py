"""Strict area and phase margins for three-strip constructions with tangent caps."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

from rational_trig import PI, add, mul, neg, scale, trig

sigma = F("1.5494009336023418793")
tolerance = F(1, 1000)  # radians for each normal, relative to fixed alpha_0
cap_slack = F(1, 100000)

def at_degree(value, function):
    return trig(scale(PI, (F(value)%360)/180), function)

def divide_nonzero_positive(a, b):
    assert b[0]>0
    return mul(a,(1/b[1],1/b[0]))

def print_fraction(x):
    return [str(y) for y in x] if isinstance(x, tuple) else str(x)

# Uniform perturbation of two cap residue phases.
assert trig(add(scale(PI,F(2,3)),(-tolerance,tolerance)),"sin")[0] > F(1,2)
sine_2 = at_degree(2,"sin")[0]
sine_8 = at_degree(8,"sin")[0]
center_gap = 90 - 3*sigma
assert 84<center_gap<90
norms = {7:scale(at_degree(3*sigma,"cos"),F(2)),
         13:scale(at_degree(6*sigma,"sin"),F(2))}
harmonics = {}
for k in (7,13):
    remainder = F(k*k+2*k-1)*tolerance / (F(3,4)*(k-1))
    residue_error = F(3*k-1)*tolerance
    residue_min=norms[k][0]-residue_error
    assert remainder < sine_2
    assert residue_error < norms[k][0]*sine_8
    excess=F(1,4*(k*k-1))*residue_min*sine_2
    assert excess > cap_slack
    harmonics[str(k)] = {
        "relative_S_minus_linear_upper":str(remainder),
        "S_angle_error_strictly_below_degrees":2,
        "exact_base_residue_norm_lower":str(norms[k][0]),
        "residue_perturbation_absolute_upper":str(residue_error),
        "residue_angle_error_strictly_below_degrees":8,
        "robust_violating_support_excess_lower":str(excess),
    }

# Independent one-axis obstruction holds over a much wider degree interval.
interval_alpha=mul(PI,(F(59,180),F(61,180)))
sinalpha=trig(interval_alpha,"sin")
sin5alpha=trig(scale(interval_alpha,F(5)),"sin")
monotone=add(scale(sin5alpha,F(-5)),sinalpha)
ratio=divide_nonzero_positive(sin5alpha,sinalpha)
assert monotone[0]>0
def one_axis_coefficient(beta):
    return add(at_degree(5*beta,"sin"),neg(mul(at_degree(beta,"sin"),ratio)))
one_axis_C=one_axis_coefficient(270+sigma)
one_axis_E=one_axis_coefficient(150+sigma)
assert one_axis_C[1]<0 and one_axis_E[0]>0

# One-cap tangent polygon.  The one-variable objective
# x/(2*(1-x^2))+x/2+(1-x)/(1+x) has denominator-positive derivative
# sign equal to P=x^4-5*x^2+8*x-2.
def polynomial(x):
    return x**4 - 5*x*x + 8*x - 2

def area(x):
    return x/(2*(1-x*x))+x/2+(1-x)/(1+x)

x0,x1 = F("0.30826229099909486615"),F("0.30826229099909486616")
assert polynomial(x0)<0<polynomial(x1)
first_component=divide_nonzero_positive((x0,x1),(2*(1-x1*x1),2*(1-x0*x0)))
third_component=divide_nonzero_positive((1-x1,1-x0),(1+x0,1+x1))
minimum_area=add(add(first_component,(x0/2,x1/2)),third_component)
assert F("0.8531920598769746041")<minimum_area[0]<=minimum_area[1]<F("0.8531920598769746045")

chosen_x=F("0.308262291")
chosen_area=area(chosen_x)
assert chosen_area<F("0.853192059877")
regular=(F(2,1),F(2,1))  # The exact regular 1-cap area is sqrt(3)/2 - cap.
# With x=tan(15°)=2-sqrt(3), bound it directly by the same rational formula.
sqrt3=(F("1.73205080756887729352"),F("1.73205080756887729354"))
assert sqrt3[0]**2<3<sqrt3[1]**2
reg_x=add(regular,neg(sqrt3))
reg_a=add(add(divide_nonzero_positive(reg_x,(2*(1-reg_x[1]**2),2*(1-reg_x[0]**2))),
              scale(reg_x,F(1,2))),
          divide_nonzero_positive((1-reg_x[1],1-reg_x[0]),(1+reg_x[0],1+reg_x[1])))
assert reg_a[0]-chosen_area>F("0.0024")

print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "trig_source_sha256":hashlib.sha256(Path(__file__).with_name("rational_trig.py").read_bytes()).hexdigest(),
    "cap_sigma_degrees_exact":str(sigma),
    "each_relative_normal_deviation_radians_at_most":str(tolerance),
    "each_cap_outward_slack_at_most":str(cap_slack),
    "two_harmonic_same_sign_window_center_separation_degrees":str(center_gap),
    "phase_audit":harmonics,
    "single_axis_5_harmonic_audit":{
        "alpha_degrees_interval_exact":["59","61"],
        "derivative_cos5alpha_minus_cosalpha_positive_lower":str(monotone[0]),
        "C_residue_coefficient":print_fraction(one_axis_C),
        "E_residue_coefficient":print_fraction(one_axis_E),
    },
    "one_cap":{
        "optimal_tangent_quarter_angle_interval_exact":[str(x0),str(x1)],
        "optimum_area_interval_exact":print_fraction(minimum_area),
        "regular_area_interval_exact":print_fraction(reg_a),
        "explicit_rational_trial_tangent_quarter_angle":str(chosen_x),
        "explicit_rational_trial_area_exact":str(chosen_area),
    },
    "passed":True,
},indent=2))
