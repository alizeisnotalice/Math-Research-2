import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time

def transforms(n,v):
 x,y=v%n,v//n
 return [(x,y),(n-1-x,y),(x,n-1-y),(n-1-x,n-1-y),(y,x),(n-1-y,x),(y,n-1-x),(n-1-y,n-1-x)]
def canonical(n,p,q):
 return min(tuple(sorted((a[1]*n+a[0],b[1]*n+b[0]))) for a,b in zip(transforms(n,p),transforms(n,q)))
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=3);a=p.parse_args()
start=time.perf_counter();out={'binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest(),'workers':a.workers,'results':[]}
for n,A2,sq in ((4,779299174080,5272309842240),(5,7063261481362958280000,67631454448823894520000)):
 orbits=Counter(canonical(n,p,q) for p,q in itertools.combinations(range(n*n),2) if abs(p%n-q%n)+abs(p//n-q//n)>1)
 def count(pair):return pair,int(subprocess.check_output([str(a.binary),str(n),str(pair[0]),str(pair[1])],text=True))
 with ThreadPoolExecutor(max_workers=a.workers) as ex:counts=dict(ex.map(count,orbits))
 hist=Counter();distance2=0
 for (p,q),count in counts.items():
  dx,dy=p%n-q%n,p//n-q//n;hist[abs(dx)+abs(dy)]+=count*orbits[p,q];distance2+=(dx*dx+dy*dy)*count*orbits[p,q]
 assert sum(hist.values())==A2 and distance2==sq
 result={'n':n,'A2':A2,'squared_euclidean_sum':distance2,'distance_counts':dict(sorted(hist.items())),
         'pair_orbits':[{'pair':list(pair),'orbit_size':orbits[pair],'count_per_pair':count} for pair,count in counts.items()],
         'tail_at_n_over_2':str(Fraction(sum(c for d,c in hist.items() if 2*d>n),A2))}
 out['results'].append(result);print(json.dumps({k:v for k,v in result.items() if k!='pair_orbits'}),flush=True)
 a.output.write_text(json.dumps(out,indent=2)+'\n')
out['wall_seconds']=time.perf_counter()-start;a.output.write_text(json.dumps(out,indent=2)+'\n')
