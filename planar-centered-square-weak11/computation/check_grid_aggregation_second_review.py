#!/usr/bin/env python3
"""Independent dimensionless Fraction checks for the grid-aggregation second review.
Finite scalar/neighbor evidence only; no global halo enumeration.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from hashlib import sha256
import json, os, sys
assert os.path.ismount('/Volumes/data')
h=F(1);s=h/8;eps=s/4096;delta=s/2;rho=F(1,32)
source_separations=[]
for i,j in product(range(-2,3),repeat=2):
    if not(i or j):continue
    for a,b in product((-1,1),repeat=2):
        source_separations.append(max(abs(i*h+(b-a)*s/2),abs(j*h)))
assert min(source_separations)==F(7,8)
r0=min(source_separations)/2
cross=rho*(1+(h+2*delta)/(2*r0))**2
aggregate=rho*(1+h/(2*(h/2)))**2
assert cross==F(8,49)<1 and aggregate==F(1,8)<1
assert s/2+eps<h/2
# A lower sqrt2 bound proves every shrunken pair rectangle remains nonempty.
assert F(7,5)**2<2 and F(10,7)**2>2
assert (F(7,5)-1)*s>4*eps
assert F(7,5)*s>4*eps
assert 4*(2*F(10,7)-1)<8
# K*s^2=1/2, hence the universal aggregate loss and fixed deficit.
loss=24*F(1,2)/4096
assert loss==F(3,1024)
gap=F(3,2)-F(10,7)-loss
assert gap==F(1,14)-F(3,1024)>F(1,16)
levels=[]
for n in (3,7,12):
    h=F(1,2**n);s=h/8;K=32/(h*h);mass=2*s*s;d=F(1,2**(4*n))
    assert K*mass==1 and K*s*s==F(1,2)
    q=mass/d;assert q.denominator==1 and q==2**(2*n-5)
    levels.append({'n':n,'cells':str(K),'cell_mass':str(mass),'quantum_multiple':str(q)})
p=Path(__file__).resolve()
result={'status':'PASS','python':sys.version,'script_sha256':sha256(p.read_bytes()).hexdigest(),
 'arithmetic':'Fraction exact; no floats','seed':None,'tolerance':None,
 'normalized_neighbor_pairs':len(source_separations),'minimum_source_separation':str(min(source_separations)),
 'cross_cell_tail_at_r0':str(cross),'aggregate_tail_at_h_over_2':str(aggregate),
 'total_inner_rectangle_loss_bound':str(loss),'deficit_rational_lower_bound':str(gap),
 'levels':levels,'scope':'finite local/scalar checks; all-scale and all-n proof separately reviewed'}
out=p.with_suffix('.json');out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
