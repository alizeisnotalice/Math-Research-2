import sys,json
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
ctx.prec=160
C=Path(sys.argv[1]).resolve();T=Path(sys.argv[2]);sys.path.insert(0,str(C))
from certify_polygon_083731 import lower_radius
from select_tomar_083731_scores import radius_rows
q=lambda x:arb(x.numerator)/x.denominator
pi=arb.pi();s=F('0.999999006');g=F(120696470013529293,2**57);root=json.loads((T/'part0000.json').read_text());lo=list(map(F,root['lo']));hi=list(map(F,root['hi']));rad=radius_rows(json.loads((C/'five-body-radial-affine.json').read_text()),F('.83750'));out=[]
for i,row in enumerate(rad):
 k=3*i;edge=min(-lo[k],-lo[k+1],hi[k],hi[k+1]);a=q(s/2);b=q(s)*lower_radius(row['n'],pi);value=pi*(a*a+b*b)/2+q(edge)*(a+b);assert value>q(g); assert F(row['integer_radius'],2**44)<edge and row['integer_area_lower']>=120696470013529293;out.append(dict(n=row['n'],root_edge=str(edge),scaled_exterior_lower=str(value.lower()),radius=row['integer_radius'],radial_area_lower=row['integer_area_lower']))
print(json.dumps(out,indent=2))
