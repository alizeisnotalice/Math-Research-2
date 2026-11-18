"""Verify on the R214 geometry that the entire reflected normal band avoids five rounded cells."""
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import sys
import platform
import flint
from flint import arb, arb_series, ctx

parser = argparse.ArgumentParser()
parser.add_argument('computation', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
comp = args.computation.resolve()
src = comp/'b1-tangent-curvature'
sys.path.insert(0, str(src))
from certify_b1_classification import Geometry, PINS, POSITIVE, SELECT, counterexample, dot, perp, rat, rotate, unit

G = Geometry(comp)
ctx.prec = 192
ctx.cap = 3
pi = arb.pi()
fractional = (69,626,1209,1808,2391)
lo,hi = Fraction(4192,10000),Fraction(4383,10000)

# Handle reflection and the antipodal width identity before testing a normal.
def half_circle_bin(u):
    angle=arb.atan2(u[1],u[0]); value=angle*2880/pi
    if value<0:
        value+=2880
    assert 0<value<2880, (angle,value,u)
    for k in fractional:
        assert value<k or value>k+1, (value,k)
    return value

ah=half_circle_bin(G.normal)
sigma=rat(Fraction(15497,10000))*pi/180
caps=[]
for k in range(0,360,60):
    for sign in (1,-1):
        if (k,sign) in ((0,1),(240,-1)):
            continue
        choice = 3 if SELECT[len(caps)]==0 else 5
        u=unit(pi/6+sigma+choice*pi/3)
        caps.append(str(half_circle_bin(rotate((u[0],sign*u[1]),-k*pi/180))))

def normal_band(a,b):
    w=arb_series([rat(a).union(rat(b)),1,0]); x,y=G.ns['points'](w)[0][2]
    q=G.qa((x,y)); v=(q[0][1],q[1][1]); speed=dot(v,v).sqrt()
    assert speed>0
    return half_circle_bin(G.qa(tuple(-z/speed for z in perp(v))))

band=None
N=8192
for j in range(N):
    val=normal_band(lo+(hi-lo)*j/N,lo+(hi-lo)*(j+1)/N)
    band=val if band is None else band.union(val)

assert 75 < band < 85
expected_base = '584e55bf2b54eaf3273f1d77c03640a8294ae3af9b1b3e1a08daf75f36cbb4dd'
expected_extension = 'cebb1ca862eb6484625da2f053946299604757a6371ce8955dc778291b0e051f'
assert hashlib.sha256((src/'certify_b1_classification.py').read_bytes()).hexdigest() == expected_base
assert hashlib.sha256((src/'extension.json').read_bytes()).hexdigest() == expected_extension
record=counterexample(G,56,POSITIVE[56]); coeff={k:arb(v) for k,v in record['rho_solution'].items()}; coeff['626']=arb(POSITIVE[56][1]); endpoints={}
for k in fractional:
    r=coeff[str(k)]; fraction=r if k//960!=1 else 1-r
    midpoint=arb(2*k+1)/2; half=(fraction*(pi/5760).sin()).asin()*2880/pi
    assert 0<half<arb(1)/2
    endpoints[str(k)]={'centered_interval':'empty' if k//960==1 else 'filled','left_bin':str(midpoint-half),'right_bin':str(midpoint+half)}

out={'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':platform.python_version(),'python_flint':flint.__version__,'precision_bits':ctx.prec,'upstream_base_sha256':hashlib.sha256((src/'certify_b1_classification.py').read_bytes()).hexdigest(),'upstream_extension_sha256':hashlib.sha256((src/'extension.json').read_bytes()).hexdigest(),'upstream_geometry':PINS,'AH_bin':str(ah),'ten_other_bins':caps,'all_reflected_bins':str(band),'band_subintervals':N,'w_interval':[str(lo),str(hi)],'replayed_witness':record,'half_circle_replacements':endpoints,'positive_curvature_components_on_circle':15}
args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('AH',ah,'QA',band,'arcs',15)
