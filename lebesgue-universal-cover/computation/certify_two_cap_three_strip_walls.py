"""Rigorous area and two-frequency checks for arbitrary two-cap strip designs."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import hashlib
import json
import sys

source=Path(__file__).with_name("certify_gibbs_area_refined.py")
raw=source.read_bytes()
parts=raw.decode().split("width = F(1,cells)")
assert len(parts)==2
sys.argv=[str(source),"1.5494009336023418793","2","1"]
context={"__file__":str(source)}
exec(compile(parts[0],str(source),"exec"),context)
I=context["I"]
trig=context["trig"]
pi=context["pi_machin"]

def unit_degree(a):
    a=I(a)
    radians=trig.mul(pi,(a.lo/180,a.hi/180))
    return tuple(I(*trig.trig(radians,k)) for k in ("cos","sin"))

def tan(a):
    c,s=unit_degree(a)
    return s/c

def sec_squared(a):
    c,_=unit_degree(a)
    return 1/(c*c)

def csc_squared(a):
    _,s=unit_degree(a)
    return 1/(s*s)

def area_same(g):
    return (tan(g/2)+3*tan(g/6)+4*tan((180-g)/4))/4

def same_derivative(g):
    return (sec_squared(I(g)/2)+sec_squared(I(g)/6)-2*sec_squared((180-I(g))/4))/8

def area_distinct(g):
    return tan(g/2)/2+tan(g/4)+1/(2*tan(g))

def distinct_derivative(g):
    return (sec_squared(I(g)/2)+sec_squared(I(g)/4)-2*csc_squared(I(g)))/4

def bracket(label,area,der,lo,hi,enclosure_lower,enclosure_upper):
    lr,hr=F(lo),F(hi)
    negative,positive=der(lr),der(hr)
    assert negative.hi<0<positive.lo
    interval=area(I(lr,hr))
    assert F(enclosure_lower)<interval.lo<=interval.hi<F(enclosure_upper)
    return {"root_interval_degrees_exact":[str(lr),str(hr)],"negative_derivative":negative.out(),
            "positive_derivative":positive.out(),"area_interval_exact":interval.out()}

same=bracket("same",area_same,same_derivative,"70.29138887277320917565","70.29138887277320917567",
             "0.8504785984619519762","0.8504785984619519765")
distinct=bracket("distinct",area_distinct,distinct_derivative,"64.05138508747993758508","64.05138508747993758511",
                 "0.8430460106436789621","0.8430460106436789624")

# The symmetric bisectors are not competitive outside [60°,68°].
endpoint_bounds={}
for angle in (60,68):
    interval=area_distinct(I(angle))
    assert interval.lo>F(".8449") if angle==68 else interval.lo>F(".8452")
    endpoint_bounds[str(angle)]=interval.out()
assert distinct_derivative(60).hi<0<distinct_derivative(68).lo
current_upper=F(".844093590940968945")
large_cut_gap_min=(1+2*tan(I(F(45,2)))+I(pi[0],pi[1])/2)/4
large_uncut_gap_min=area_distinct(I(45))
assert large_cut_gap_min.lo>current_upper and large_uncut_gap_min.lo>current_upper

# Counterexample on every mirror-bisector config with beta in [30,34]
# except beta=30. Uniform robust inequality for a phase t in [52,68],
# evaluated in 16 x 128 cells. Odd multiples of beta supply such a t.
B=F(3,2)
num_beta,num_t=16,128
min_margin=F(1000)
smallest_A1=F(1000)
largest_q=F(0)
for i in range(num_beta):
    beta=I(30+F(4*i,num_beta),30+F(4*(i+1),num_beta))
    c,_=unit_degree(beta)
    for j in range(num_t):
        t=I(52+F(16*j,num_t),52+F(16*(j+1),num_t))
        ct,st=unit_degree(t)
        c3,s3=unit_degree(3*t)
        A1,A3=ct+c,c3+c
        assert A1.lo>0>A3.hi and (ct-c).hi<0 and (c3-c).hi<0
        q=-B*A3/A1
        max_sin3=max(abs(s3.lo),abs(s3.hi))
        max_sin1=max(abs(st.lo),abs(st.hi))
        gap=B*(1-A3/A1)*(1-4*q*q)-I(max_sin1)/c-3*B*q*I(max_sin3)/(c*(1-q*q).sqrt())
        assert gap.lo>F(1,10), (i,j,gap.out())
        min_margin=min(min_margin,gap.lo)
        smallest_A1=min(smallest_A1,A1.lo)
        largest_q=max(largest_q,q.hi)

# beta/pi = p/q in lowest terms. Odd multiples traverse q residues on
# a translated uniform grid. If q>=23, its step is <16 degrees.
assert F(360,23)<16
small_denominators={}
for denominator in range(1,23):
    for numerator in range(1,denominator):
        frac=F(numerator,denominator)
        if gcd(numerator,denominator)!=1 or not (F(1,6)<=frac<=F(17,90)):
            continue
        choices=[]
        for j in range(denominator):
            d=2*j+1
            phase=F(180)*frac*d%360
            unsigned=min(phase,360-phase)
            if 52<=unsigned<=68:
                choices.append((d,phase))
        if frac==F(1,6):
            assert not choices
            small_denominators[str(frac)]={"regular_exception":True}
        else:
            assert choices and choices[0][0]>=3
            small_denominators[str(frac)]={"frequency_d":choices[0][0],"phase_degrees_exact":str(choices[0][1])}

print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "coordinate_source_sha256":hashlib.sha256(raw).hexdigest(),
    "both_caps_same_gap_geometric_minimum":same,
    "caps_different_gaps_geometric_minimum":distinct,
    "large_cut_gap_area_lower":str(large_cut_gap_min.lo),
    "large_uncut_gap_area_lower":str(large_uncut_gap_min.lo),
    "bisector_family_endpoint_areas":endpoint_bounds,
    "mirror_bisector_two_frequency":{
        "beta_degrees_interval":[30,34],"phase_degrees_interval":[52,68],
        "interval_boxes":num_beta*num_t,
        "least_normalized_support_margin_lower":str(min_margin),
        "smallest_A1_lower":str(smallest_A1),
        "largest_bound_on_abs_cos_psi":str(largest_q),
        "small_denominators_explicit_cases":small_denominators,
    },
    "passed":True,
},indent=2))
