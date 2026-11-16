from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,platform,time

def ceilq(x):return -((-x.numerator)//x.denominator)
def floorq(x):return x.numerator//x.denominator
def cnt(lo,hi):return max(0,hi-lo+1)
def main(inp,out):
 start=time.perf_counter();config=json.loads(inp.read_text());records=[];counts=dict(cap_pairs=0,forward_differences=0,interval_overlap_checks=0)
 for n in config['n_values']:
  e=F(1,n*n);h=F(n+1,n*n);w=1-e;t=F(n,n+1);grid=[i*h for i in range(n)];ak=lambda k:(k-1)*h+e
  ff=lambda k,l:F(k*l,n*n)-t*ak(k)*ak(l)
  for k in range(1,n+1):
   for l in range(1,n+1):
    assert ff(k,l)<=1-t;counts['cap_pairs']+=1
    if k<n:assert ff(k+1,l)-ff(k,l)==(1-F(l,n))/(n*n);counts['forward_differences']+=1
  for R in [F(1),1+e,F(2)]:
   events=sorted(set(grid+[p-R for p in grid]+[-e,F(1),-e-R,1-R]));samples=events+[(a+b)/2 for a,b in zip(events,events[1:])]
   for left in samples:
    k=sum(left<=p<=left+R for p in grid)
    if k:
     overlap=max(F(0),min(left+R,F(1))-max(left,-e));assert overlap>=ak(k);counts['interval_overlap_checks']+=1
  a=F(n+2,n*n);k=sum(-a<=p<=1-a for p in grid);assert k==n-1;assert a-e==h;assert F(k,n)+t*h==1
  U=(1+e)**2+F(3-6*e,n+1)
  # Exact Stieltjes layer-cake integral for left-unit-strip cumulative caps.
  # Integrating 1-k_n(a)/n over 0<=a<=1 gives w/2.
  integral=1-sum(1-p for p in grid)/n;assert integral==w/2
  records.append(dict(n=n,density=str(t),sharp_unit_shift=str(a),atomization_limit_upper=str(U),left_strip_integrated_cap=str(integral)))
 cells=[]
 for n in config['atomization_n_values']:
  e=F(1,n*n);w=1-e;t=F(n,n+1);lo=F(-1,2);hi=w+F(1,2);jlo=-e;jhi=F(1)
  for m in config['mesh_denominators']:
   delta=F(1,m);il=ceilq(lo/delta);iu=floorq(hi/delta-1);nc=cnt(il,iu)
   bl=max(il,ceilq(jlo/delta-1));bu=min(iu,floorq(jhi/delta));nb=cnt(bl,bu)
   outside_cells=nc*nc-nb*nb;cell_area=outside_cells*delta*delta;target=(2-e)**2-(1+e)**2
   assert 0<=cell_area<=target
   cells.append(dict(n=n,mesh=m,complete_outside_cells=outside_cells,complete_cell_area=str(cell_area),residual_singleton_halo_area_lower=str(t*cell_area),limiting_residual_lower=str(t*target),missing_boundary_cell_area=str(target-cell_area)))
 result=dict(status='PASS_EXACT_BOREL_SCALARS_AND_ATOMIC_CELL_COUNTS',python=platform.python_version(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),counts=counts,records=records,atomic_cell_counts=cells,elapsed_seconds=time.perf_counter()-start,scope='Finite rational checks of all-n analytic formulas and exact singleton-halo area lower counts for prescribed canonical grids. No full halo computation, no actual maximum eligible mass evaluation for those large atomizations, no finite-atomic3.5 claim.')
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.input,a.output)
