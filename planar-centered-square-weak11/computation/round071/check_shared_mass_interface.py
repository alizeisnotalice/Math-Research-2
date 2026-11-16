from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse,hashlib,json,platform,time

def frac(v):return F(str(v))
def submasks(n):return range(1,1<<n)
def ids(mask,n):return [i for i in range(n) if mask>>i&1]
def mass(mask,m):return sum((m[i] for i in ids(mask,len(m))),F(0))
def bbox(mask,p):
 s=[p[i] for i in ids(mask,len(p))]
 return min(x for x,y in s),max(x for x,y in s),min(y for x,y in s),max(y for x,y in s)
def span(b):return max(b[1]-b[0],b[3]-b[2])
def closed_rect(q,b):return b[0]<=q[0]<=b[1] and b[2]<=q[1]<=b[3]
def dist(x,q):return max(abs(x[0]-q[0]),abs(x[1]-q[1]))
def halo_subsets(x,p,m):
 for k in submasks(len(p)):
  if mass(k,m)>4*max(dist(x,p[i]) for i in ids(k,len(p)))**2:return True
 return False
def prefix_nonhalo(x,p,m):
 return all(sum((m[i] for i in range(len(p)) if dist(x,p[i])<=r),F(0))<=4*r*r for r in set(dist(x,q) for q in p))
def menus_halo(x,p,m):
 return any(sum((m[i] for i in range(len(p)) if dist(x,p[i])<=r),F(0))>4*r*r for r in set(dist(x,q) for q in p) if r<F(1,2))

def main(inp,out):
 start=time.perf_counter();o=json.loads(inp.read_text());S=[tuple(map(frac,a)) for a in o['S_positions']];sm=list(map(frac,o['S_masses']));Q=[tuple(map(frac,a)) for a in o['Q_positions']];P=S+Q;xs=[tuple(map(frac,a)) for a in o['receivers']]
 closures={}
 for k in submasks(len(P)):
  b=bbox(k,P);cl=sum(1<<i for i,p in enumerate(P) if closed_rect(p,b));closures[cl]=span(b)**2
 counts=dict(mass_vectors=0,all_source_subsets=0,admissible_vectors=0,residual_receiver_checks=0,admissible_full_menu_checks=0,admissible_deletion_checks=0)
 for ys in product(map(frac,o['outside_mass_menu']),repeat=len(Q)):
  m=sm+list(ys);counts['mass_vectors']+=1
  raw_rows=[not(mass(k,m)>1 and span(bbox(k,P))**2<mass(k,m)) for k in submasks(len(P))]
  raw=all(raw_rows)
  caps=all(mass(k,m)<=max(F(1),d2) for k,d2 in closures.items())
  assert raw==caps;counts['all_source_subsets']+=(1<<len(P))-1
  if caps:counts['admissible_vectors']+=1
  for x in xs:
   nh=prefix_nonhalo(x,Q,list(ys));assert nh==not_halo(x,Q,list(ys));counts['residual_receiver_checks']+=1
   if caps:
    full=halo_subsets(x,P,m);menu=menus_halo(x,P,m);assert full==menu;counts['admissible_full_menu_checks']+=1
    assert (full and not halo_subsets(x,Q,list(ys)))==(menu and nh);counts['admissible_deletion_checks']+=1
 # Strict threshold and zero-distance/tie edge cases.
 assert prefix_nonhalo((F(1,2),F(0)),[(F(0),F(0))],[F(1)])
 assert not prefix_nonhalo((F(0),F(0)),[(F(0),F(0))],[F(1)])
 assert prefix_nonhalo((F(0),F(0)),[(F(0),F(0))],[F(0)])
 assert not prefix_nonhalo((F(0),F(0)),[(F(-1,4),F(0)),(F(1,4),F(0))],[F(1,4),F(1,4)])
 # Matrix-only weak-duality diagnostics; these are not geometric examples.
 alpha=F(1,2);lam=[F(1,2),F(1,2)];d=[F(1,2),F(1,2)]
 assert all(alpha>=l for l in lam)
 bound=alpha-sum(l*v for l,v in zip(lam,d));assert bound==0
 for y in [(F(3,4),F(0)),(F(0),F(3,4))]:assert sum(y)<=1
 result=dict(status='PASS_EXACT_FINITE_INTERFACE_CHECKS',python=platform.python_version(),input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),counts=counts,distinct_bbox_closures=len(closures),strict_edge_fixtures=4,matrix_only_dual_upper_bound=str(bound),elapsed_seconds=time.perf_counter()-start,scope='Fixed seven-position exhaustive algebra checks, including zero external mass and forced-zero inside-bbox column; exact residual complete-prefix and finite witness-menu equivalence. Matrix diagnostic is not an actual geometric joint obstruction. No LP optimization, halo area or C2 bound.')
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
def not_halo(x,p,m):return not halo_subsets(x,p,m)
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--input',type=Path,required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();main(args.input,args.output)
