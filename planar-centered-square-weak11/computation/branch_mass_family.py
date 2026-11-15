#!/usr/bin/env python3
"""Exact branch-and-bound of an affine two-parameter periodic-mass family."""
import argparse,hashlib,heapq,json,os,platform,sys,time
from fractions import Fraction as F
from pathlib import Path
from mass_family_envelope import MassFamilyEnvelope
from periodic_halo_dominated import calculate_dominated
from periodic_halo import union_area

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--output',required=True);p.add_argument('--gap',default='1/100000000');p.add_argument('--max-leaves',type=int,default=5000);p.add_argument('--seconds',type=float,default=3600);p.add_argument('--initial',nargs=2);p.add_argument('--bits',type=int,default=40);p.add_argument('--radius',default='5');p.add_argument('--allow-zero-corners',action='store_true',help='valid upper bounds only; no strict-positivity convergence budget');args=p.parse_args()
 gap=F(args.gap)
 if gap<=0 or args.max_leaves<1 or args.seconds<=0:p.error('positive gap, max leaves and seconds required')
 start=time.monotonic();family=MassFamilyEnvelope(json.loads(Path(args.input).read_text()),args.bits,F(args.radius));init=time.monotonic()-start
 if not family.meta['all_corner_labels_strictly_positive'] and not args.allow_zero_corners:raise ValueError('strictly positive corner labels required unless explicitly overridden')
 samples={};best=None;bestpoint=None
 def sample(a,b):
  nonlocal best,bestpoint
  key=(a,b)
  if key not in samples:
   data=family.point(a,b);r=calculate_dominated(data,family.base,family.copy_factors,family.delta,bits=family.bits,R=family.R)
   lower,upper=F(r['inner_ratio']),F(r['outer_ratio']);samples[key]=(lower,upper)
   if best is None or lower>best:best,bestpoint=lower,key
 def bound(cell):return F(family.upper(*cell)['upper_ratio'])
 a0,a1,b0,b1=family.a0,family.a1,family.b0,family.b1
 point=tuple(map(F,args.initial)) if args.initial else ((a0+a1)/2,(b0+b1)/2)
 sample(*point);u=bound((a0,a1,b0,b1));assert best<=u
 heap=[(-u,0,(a0,a1,b0,b1))];finished=[];serial=divisions=0;reason='time_limit'
 while heap:
  if -heap[0][0]<best:
   while heap:
    nu,_,cell=heapq.heappop(heap);finished.append((cell,-nu))
   reason='target_gap';break
  if -heap[0][0]-best<=gap:reason='target_gap';break
  if len(heap)+len(finished)>=args.max_leaves:reason='leaf_limit';break
  if time.monotonic()-start>=args.seconds:break
  nu,_,(al,ah,bl,bh)=heapq.heappop(heap);am,bm=(al+ah)/2,(bl+bh)/2;sample(am,bm);assert samples[(am,bm)][0]<=-nu
  if (ah-al)/(a1-a0)>=(bh-bl)/(b1-b0):children=((al,am,bl,bh),(am,ah,bl,bh))
  else:children=((al,ah,bl,bm),(al,ah,bm,bh))
  for cell in children:
   z=bound(cell);serial+=1
   if z<best:finished.append((cell,z))
   else:heapq.heappush(heap,(-z,serial,cell))
  divisions+=1
 leaves=[{'cell':list(map(str,c)),'upper':str(u),'status':'pruned'} for c,u in finished]+[{'cell':list(map(str,c)),'upper':str(-nu),'status':'active'} for nu,_,c in heap]
 rectangles=[tuple(map(F,v['cell'])) for v in leaves];domain_area=(a1-a0)*(b1-b0)
 assert all(a0<=al<ah<=a1 and b0<=bl<bh<=b1 for al,ah,bl,bh in rectangles)
 assert sum(((ah-al)*(bh-bl) for al,ah,bl,bh in rectangles),F(0))==union_area(rectangles)==domain_area
 global_upper=max([best]+[F(v['upper']) for v in leaves]);assert all(F(v['upper'])<best for v in leaves if v['status']=='pruned')
 if reason=='target_gap':assert global_upper-best<=gap
 result={'input_sha256':sha(args.input),'driver_sha256':sha(__file__),'envelope_sha256':sha(Path(__file__).with_name('mass_family_envelope.py')),'lower_kernel_sha256':sha(Path(__file__).with_name('periodic_halo_dominated.py')),'python':sys.version,'platform':platform.platform(),'family':family.meta,'target_gap':str(gap),'max_leaves':args.max_leaves,'status':reason,'partition_area_and_union_verified':True,'initialization_seconds':init,'elapsed_seconds':time.monotonic()-start,'divisions':divisions,'best_point':list(map(str,bestpoint)),'lower':str(best),'upper':str(global_upper),'gap':str(global_upper-best),'samples':[{'point':list(map(str,k)),'lower':str(v[0]),'upper':str(v[1])} for k,v in sorted(samples.items())],'leaves':leaves}
 temp=Path(args.output+'.tmp');temp.write_text(json.dumps(result,indent=2)+'\n');os.replace(temp,args.output)
 print(json.dumps({k:v for k,v in result.items() if k not in ('samples','leaves')},indent=2))
if __name__=='__main__':main()
