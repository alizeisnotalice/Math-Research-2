"""Rigorous nominal area tails for the complete ordering range -52 < d < 32."""
from fractions import Fraction as F
from hashlib import sha256
import json,sys
from pathlib import Path
from flint import arb,ctx
ctx.prec=160
def q(x):
    x=F(x);return arb(x.numerator)/x.denominator
def deg(z):return q(z)*arb.pi()/180
def bounds(x):
    ans=[]
    for y in [x.lower(),x.upper()]:
        m,e=y.man_exp();ans.append(str(F(int(m))*F(2)**int(e)))
    return ans
def geometry(v):
    d=deg(v)
    gs=[deg(z) for z in (64,32,32,52,32,32,64,52)]
    slope=(0,0,-1,1,0,0,1,-1)
    gs=[x+s*d for x,s in zip(gs,slope)]
    A=sum((g/2).tan() for g in gs)/4
    deriv=sum(s/(g/2).cos()**2 for s,g in zip(slope,gs))/8
    return A,deriv
a,da=geometry(-18);b,db=geometry(2)
assert da<0<db
# Every nominal adjacent gap in the natural ordering stays below 104 degrees.
# Sorted perturbed angles are within epsilon coordinatewise of the sorted nominal
# angles, by the order statistics Lipschitz inequality, even if a pair crosses.
eps=deg('.001')
max_gap=deg(104)+2*eps
lip=((max_gap/2).cos()**(-2)-1)/8
wall=min(a.lower(),b.lower())-8*eps*lip
assert wall>q('.84454')
sample=geometry(1)[0]
assert q('.84376436370')<sample<q('.84376436371')
res={'sample_one_degree':bounds(sample),'nominal_sorted_range_degrees':['-52','32'],'inner_degrees':['-18','2'],
    'normal_perturbation_deg':'.001','endpoint_areas':list(map(bounds,(a,b))),
    'endpoint_derivatives':list(map(bounds,(da,db))),'tail_floor':bounds(wall),
    'prec':ctx.prec,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
Path(sys.argv[1]).write_text(json.dumps(res,indent=1)+'\n');print(float(wall))
