"""Exhaustive n=2,3 audit of the long-geodesic high-beta family injection."""
from itertools import permutations
from pathlib import Path
from fractions import Fraction
import json
B=Path(__file__).parent/'results'
old=json.loads((B/'geodesic-r001.json').read_text())['results']
rows=[]
for n,A2 in ((2,8),(3,143808)):
 N=n*n;d=2*n-1;adj=[]
 for v in range(N):
  x,y=v%n,v//n;adj.append({w for w,ok in ((v-1,x>0),(v+1,x<n-1),(v-n,y>0),(v+n,y<n-1)) if ok})
 family=set();images={};single_k=set()
 for u in (0,n-1,N-n,N-1):
  v=N-1-u;rest=[x for x in range(N) if x not in (u,v)]
  for tail in permutations(rest):
   order=(u,)+tail;prefix=order[:d-1]
   if v not in adj[prefix[-1]] or any(y not in adj[x] for x,y in zip(prefix,prefix[1:])):continue
   active={u};valid=True
   for x in order[1:]:
    if not adj[x]&active:valid=False;break
    active.add(x)
   if not valid:continue
   image=order[:d-1]+(v,)+order[d-1:]
   active={image[0]}
   for z in image[1:]:
    assert adj[z]&active
    active.add(z)
   for k in range(1,d-1):
    out=order[:k]+(v,)+order[k:];labels={x:N-i for i,x in enumerate(out)}
    peaks=[x for x in range(N) if all(labels[x]>labels[y] for y in adj[x])]
    assert set(peaks)=={u,v}
    M=max(labels[x] for x in adj[v]);new=[labels[x] for x in adj[v] if all(y==v or labels[x]>labels[y] for y in adj[x])]
    b=max(new,default=0);assert M==N-d+1 and b==0
    assert Fraction(M-b,labels[v]-1)==Fraction(N-2*n+2,N-k-1)
    recovered=image[:k]+(image[d-1],)+image[k:d-1]+image[d:]
    assert recovered==out and (image,k) not in images
    images[image,k]=out;family.add(out)
    if k==1:single_k.add(out)
 assert len(images)==len(family)
 cert=next(x for x in old if x['n']==n);gd=next(x for x in cert['rows'] if x['d']==d);g=int(gd['numerator'])
 assert len(family)<=(d-2)*g and len(single_k)<=g
 rows.append({'n':n,'family_count':len(family),'fixed_k1_count':len(single_k),'distinct_single_peak_images':len({x for x,k in images}),'all_geodesic_prefix_count':g,'max_inverse_load':d-2,'family_probability':str(Fraction(len(family),A2))})
print(json.dumps({'verified':True,'rows':rows},indent=2))
