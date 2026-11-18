"""Compare a freely placed seven-body hull with the hexagon-constrained optimum."""
from fractions import Fraction as Q
import json
from pathlib import Path
import platform
import flint
from flint import arb, ctx
from certify_nine import (ar, append_point, arc_body, fourier_body, hull,
                          arithmetic, FOURIER_TERMS, PI_QUALITY)

ctx.prec = PI_QUALITY
base = Path(__file__).parent
k = 2048
source = json.loads((base/"056-free-seven.json").read_bytes())
params = [ar(s) for s in source["active_parameters"]]
assert len(params) == 17
pi = arb.pi()
assert pi < ar(Q(22,7))
points = []
for i in range(8*k):
    angle = 2*pi*i/(8*k)
    sn,cs = angle.sin_cos()
    append_point(points,cs/2,sn/2)
arc_body(points,k,pi,3,arb(0),*params[:2])
for n,j in ((5,2),(7,5),(9,8)):
    arc_body(points,k,pi,n,*params[j:j+3])
for i, subdivisions in enumerate((8*k,16*k)):
    j = 11+3*i
    fourier_body(points,subdivisions,pi,FOURIER_TERMS[i],*params[j:j+3])
result = arithmetic(hull(points),k)
hex_lower = Q(json.loads((base/"053-hex-fourier.json").read_bytes())["lower"])
gap = hex_lower-Q(result["upper"])
assert gap > Q(".00013672")
assert Q(result["upper"]) < Q(".83796232")
assert hex_lower/(1+90*Q(".000000905"))**2 > Q(result["upper"])
print(json.dumps({"python":platform.python_version(),"python_flint":flint.__version__,
                  "precision":ctx.prec,"arc_subdivisions":k,"sample_count":len(points),
                  "area":result,"hex_lower_minus_free_upper":str(gap)},indent=2))
