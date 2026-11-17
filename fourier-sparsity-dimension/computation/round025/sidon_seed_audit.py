from pathlib import Path
from itertools import combinations
import json,time,hashlib,platform
start=time.monotonic();MOD={2:7,3:11,4:19,5:37,6:67}
def mul(a,b,k):
 o=0
 while b:
  if b&1:o^=a
  b>>=1;a<<=1
  if a&(1<<k):a^=MOD[k]
 return o
def trace(a,k):
 t=0
 for _ in range(k):t^=a;a=mul(a,a,k)
 assert t in (0,1)
 return t
def basis(vs):
 p={};out=[]
 for a in vs:
  orig=a
  while a:
   b=a.bit_length()-1
   if b in p:a^=p[b]
   else:p[b]=a;out.append(orig);break
 return out
def span(bs):
 s=[0]
 for b in bs:s+=[x^b for x in s]
 return s
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
def bent(m,N,q):return all(abs(x)==q for x in fwht([1-2*((m>>i)&1) for i in range(N)]))
rows=[]
for k in range(2,7):
 q=1<<k;cubes=[mul(mul(w,w,k),w,k) for w in range(q)]
 B=basis([(1<<q)-1]+[sum(trace(mul(1<<j,w,k),k)<<w for w in range(q)) for j in range(k)]+[sum(trace(mul(1<<j,cubes[w],k),k)<<w for w in range(q)) for j in range(k)])
 labels=[sum(((b>>w)&1)<<j for j,b in enumerate(B)) for w in range(q)]
 differences=[labels[u]^labels[v] for u,v in combinations(range(q),2)]
 assert len(set(differences))==len(differences) and 0 not in differences
 # Direct derivative check of the cubing graph, independent of trace code.
 for s in range(1,q):
  values={}
  for w in range(q):values.setdefault(cubes[w]^cubes[w^s],[]).append(w)
  assert all(len(v)==2 and v[0]^v[1]==s for v in values.values())
 row={'k':k,'q':q,'field_modulus':MOD[k],'seed_dimension':len(B),'seed_basis_hex':[hex(x) for x in B],'unordered_pairs':len(differences),'distinct_nonzero_pair_characters':True}
 if k==2:
  N=q*q;lift=lambda h:sum(((h>>(x>>k))&1)<<x for x in range(N));H=[lift(h) for h in span(B)]
  candidates=[m for m in range(1<<N) if bent(m,N,q)]
  accepted=[m for m in candidates if all(bent(m^h,N,q) for h in H)]
  def is_mm(m):
   slopes=[]
   for w in range(q):
    vals=[1-2*((m>>(z+(w<<k)))&1) for z in range(q)];W=fwht(vals);sup=[a for a,t in enumerate(W) if t]
    if len(sup)!=1 or abs(W[sup[0]])!=q:return False
    slopes+=sup
   return len(set(slopes))==q
  assert len(candidates)==896 and len(accepted)==384 and all(is_mm(m) for m in accepted)
  g=sum((((z&w).bit_count()%2)^((z&1)*((z>>1)&1)))<<(z+(w<<k)) for z in range(q) for w in range(q))
  affineW=span([15,10,12]);assert all(bent(g^lift(h),N,q) for h in affineW) and not is_mm(g)
  witness=next(h for h in span(B) if not bent(g^lift(h),N,q))
  row.update({'all_truth_tables':1<<N,'bent_count':len(candidates),'seed_compatible_count':len(accepted),'compatible_masks_hex':[hex(x) for x in accepted],'all_compatible_are_MM':True,'affine_seed_counterexample_hex':hex(g),'detecting_seed_phase_hex':hex(witness)})
 rows.append(row)
out={'rows':rows,'arithmetic':'exact GF2 fields and integer Walsh','python':platform.python_version(),'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print([{k:v for k,v in r.items() if k!='seed_basis_hex'} for r in rows])
