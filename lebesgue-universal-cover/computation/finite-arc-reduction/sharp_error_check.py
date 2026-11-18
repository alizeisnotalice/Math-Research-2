"""Reproduce the sector extremizer's moment/reachability identities and range enclosures."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform

import flint
from flint import arb, ctx
import sympy as sp

parser = argparse.ArgumentParser()
parser.add_argument('output',type=Path)
args = parser.parse_args()
ctx.prec = 192

# The center is 0, the sector [-x,x], and sin(ell)=sin(x)/2.
x,ell=sp.symbols('x ell',real=True)
s=sp.symbols('s',real=True)
def integrate_sign(weight):
    return sp.integrate(weight,(s,-x,x))-2*sp.integrate(weight,(s,-ell,ell))
assert sp.simplify(integrate_sign(sp.sin(s)))==0
assert sp.simplify(integrate_sign(sp.cos(s)).subs(sp.sin(ell),sp.sin(x)/2))==0
peak=sp.integrate(sp.sin(-s),(s,-x,0))-2*sp.integrate(sp.sin(-s),(s,-ell,0))
assert sp.simplify(peak-(2*sp.cos(ell)-1-sp.cos(x)))==0
c=sp.symbols('c',real=True)
assert sp.simplify((sp.sqrt(3+c*c)-1-c)*(sp.sqrt(3+c*c)+1+c)-(2-2*c))==0
series=sp.series(sp.sqrt(4-sp.sin(x/2)**2)-1-sp.cos(x/2),x,0,6)
assert series.coeff(x,2)==sp.Rational(1,16)

rows=[]
for n in [3,12,128,512,1024,2048,4096,8192,16384]:
    d=arb.pi()/n
    c=(d/2).cos()
    radius=(1-c)*2/((3+c*c).sqrt()+1+c)
    ell=((d/2).sin()/2).asin()
    assert (radius-(2*ell.cos()-1-c)).contains(0)
    assert 0<radius<((d/4).tan())**2
    rows.append({'half_circle_cells':n,'maximum_arcs':2*n+1,'error_enclosure':str(radius)})
args.output.write_text(json.dumps({'python':platform.python_version(),'python_flint':flint.__version__,'sympy':sp.__version__,'precision_bits':ctx.prec,'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'sharp_error_series':str(series),'grid_radius_checks':rows},indent=2)+'\n')
print('symbolic identities passed;',len(rows),'Arb enclosures')
