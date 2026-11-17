from pathlib import Path
from math import comb
from collections import Counter
import json,time,hashlib,platform,random
start=time.monotonic()
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
def rank(vs):
 p={}
 for x in vs:
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;break
 return len(p)
rows=[]
for k,ell in [(3,1),(3,2),(4,2),(4,3),(5,2),(5,3),(6,2)]:
 q=1<<k;N=q*q
 truth=lambda fn:sum(int(fn(x&(q-1),x>>k))<<x for x in range(N))
 gens=[truth(lambda z,w,u=u:(w&u)==u) for u in range(q) if u.bit_count()<=k-ell]
 for j in range(k):
  for u in range(1<<j):
   if u.bit_count()<=k-ell-1:gens.append(truth(lambda z,w,u=u,j=j:((z>>j)&1) and (w&u)==u))
 D=2*sum(comb(k,d) for d in range(k-ell+1))-1;assert len(gens)==rank(gens)==D
 target=(1<<ell)*q;base=truth(lambda z,w:(z&w).bit_count()%2)
 hist=Counter()
 if D<=21:
  word=0;hist[0]=1
  for i in range(1,1<<D):
   j=(i&-i).bit_length()-1;word^=gens[j];hist[word.bit_count()]+=1
  assert min(x for x in hist if x)>0 and min(x for x in hist if x)==target
 rng=random.Random(2026112400+10*k+ell)
 for _ in range(64):
  mask=rng.randrange(1<<D);word=0
  for j,g in enumerate(gens):
   if mask>>j&1:word^=g
  assert not word or word.bit_count()>=target
  assert all(abs(v)==q for v in fwht([1-2*(((base^word)>>x)&1) for x in range(N)]))
 witness=truth(lambda z,w:(w&((1<<(k-ell))-1))==((1<<(k-ell))-1));assert witness.bit_count()==target
 rows.append({'k':k,'ell':ell,'dimension':D,'minimum_weight_target':target,'full_weight_enumeration_count':sum(hist.values()),'weight_histogram':dict(sorted(hist.items())),'sampled_bent_phases':64,'seed':2026112400+10*k+ell,'explicit_equality_witness':True})
out={'scope':'exact weights if dimension<=21; otherwise specified samples only; all-dimensional theorem separate','rows':rows,'python':platform.python_version(),'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print([{k:v for k,v in r.items() if k!='weight_histogram'} for r in rows])
