"""All 3D kernels: coset occupancy profiles for the eight-dimensional dual block."""
import itertools, multiprocessing as mp, json, os, time
from pathlib import Path
from collections import Counter
T=[0,2,4,6,65,67,68,70,74,75,88,89,104,105,120,121];D=T+[x+128 for x in T]
def scan(worker):
 counts=Counter(); total=0
 for piv in itertools.combinations(range(8),3):
  slots=[(i,j) for i,p in enumerate(piv) for j in range(p+1,8) if j not in piv]
  for assignment in range(worker,1<<len(slots),12):
   rows=[1<<p for p in piv]
   for h,(i,j) in enumerate(slots):rows[i]|=((assignment>>h)&1)<<j
   hits=Counter()
   for x in D:
    for p,v in zip(piv,rows):
     if x>>p&1:x^=v
    hits[x]+=1
   profile=tuple(sorted(hits.values()));counts[profile]+=1;total+=1
 return os.getpid(),total,counts
if __name__=='__main__':
 start=time.monotonic()
 with mp.Pool(12) as pool:results=pool.map(scan,range(12),1)
 c=Counter();
 for pid,n,counts in results:c.update(counts)
 assert sum(c.values())==97155 and c[(1,)*32]>0
 out=dict(workers=[(p,n) for p,n,_ in results],total=sum(c.values()),injective=c[(1,)*32],profiles={str(k):v for k,v in c.items()},elapsed=time.monotonic()-start)
 Path(__file__).with_name('result.json').write_text(json.dumps(out));print(out['workers'],out['total'],out['injective'],out['elapsed'])
