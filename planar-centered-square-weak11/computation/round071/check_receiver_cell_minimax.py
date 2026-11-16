from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse,json,hashlib,platform,time

def main(inp,out):
 start=time.perf_counter();o=json.loads(inp.read_text());P=[tuple(map(F,p)) for p in o['source_positions']];cells=[tuple(map(F,c)) for c in o['cells']];radii=list(map(F,o['radii']));counts=dict(subset_cell_minima=0,exact_open_rectangle_intersections=0,cell_maxima=0)
 for mask in range(1,1<<len(P)):
  T=[p for i,p in enumerate(P) if mask>>i&1];a=min(x for x,y in T);b=max(x for x,y in T);c=min(y for x,y in T);d=max(y for x,y in T)
  for u,v,s,t in cells:
   delta=max((b-a)/2,(d-c)/2,b-v,u-a,d-t,s-c,F(0))
   # Independent finite minimization: every one-dimensional absolute-value
   # minimizer is represented among endpoints or clamped pair midpoints.
   xs={u,v};ys={s,t}
   xs.update(min(v,max(u,(p[0]+q[0])/2)) for p in T for q in T)
   ys.update(min(t,max(s,(p[1]+q[1])/2)) for p in T for q in T)
   direct=min(max(max(abs(x-px),abs(y-py)) for px,py in T) for x,y in product(xs,ys))
   assert direct==delta;counts['subset_cell_minima']+=1
   maxdist=max(max(abs(x-px),abs(y-py)) for x,y in product([u,v],[s,t]) for px,py in T)
   assert maxdist==max(abs(u-b),abs(v-a),abs(s-d),abs(t-c));counts['cell_maxima']+=1
   for r in radii:
    inter=(b-r<a+r and d-r<c+r and v>b-r and u<a+r and t>d-r and s<c+r)
    assert inter==(delta<r);counts['exact_open_rectangle_intersections']+=1
 result=dict(status='PASS_EXACT_RECEIVER_CELL_GEOMETRY',python=platform.python_version(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),counts=counts,elapsed_seconds=time.perf_counter()-start,scope='Exact rational finite minimax and open-receiver/closed-cell intersection checks, including degenerate cells and source subsets; no actual halo area or new C2 bound.')
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.input,a.output)
