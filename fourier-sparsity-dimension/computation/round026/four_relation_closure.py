from pathlib import Path
from itertools import combinations
import json,time,hashlib,platform
start=time.monotonic()
def basis(vs):
 p={};out=[]
 for orig in vs:
  x=orig
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;out.append(orig);break
 return out
def span(bs):
 a={0}
 for b in bs:a|={x^b for x in tuple(a)}
 return a
def perp(bs,q):return basis(x for x in range(1<<q) if all((x&b).bit_count()%2==0 for b in bs))
def group_pairs(H,q):
 labels=[sum(((h>>w)&1)<<j for j,h in enumerate(H)) for w in range(q)]
 groups={}
 for u,v in combinations(range(q),2):
  key=labels[u]^labels[v];assert key
  groups.setdefault(key,[]).append((u,v))
 return sorted(sorted(g) for g in groups.values())
def inspect(D,k):
 q=1<<k;words=span(D);small=[w for w in words if w.bit_count()==4];D4=basis(small);H=perp(D,q);Hb=perp(D4,q)
 affine=[(1<<q)-1]+[sum(((w>>j)&1)<<w for w in range(q)) for j in range(k)]
 assert all((a&d).bit_count()%2==0 for a in affine for d in D)
 assert group_pairs(H,q)==group_pairs(Hb,q)
 parents=list(range(q))
 def root(x):
  while parents[x]!=x:x=parents[x]
  return x
 for v in small:
  support=[w for w in range(q) if v>>w&1]
  for w in support[1:]:parents[root(w)]=root(support[0])
 comps={}
 for w in range(q):comps.setdefault(root(w),[]).append(w)
 components=list(comps.values())
 for C in components:
  mask=sum(1<<w for w in C);assert all((mask&d).bit_count()%2==0 for d in D4)
 return {'k':k,'dual_basis':D,'dual_dimension':len(D),'four_relation_rank':len(D4),'seed_dimension':len(H),'closure_dimension':len(Hb),'pair_partition_equal':True,'pair_character_classes':len(group_pairs(H,q)),'components':components}
k=3;q=8;dual=perp([(1<<q)-1]+[sum(((w>>j)&1)<<w for w in range(q)) for j in range(k)],q)
subspaces={tuple([0])};front=[{0}];whole=span(dual)
while front:
 S=front.pop()
 for v in whole-S:
  T=S|{x^v for x in S};key=tuple(sorted(T))
  if key not in subspaces:subspaces.add(key);front.append(T)
assert len(subspaces)==67
smallrows=[inspect(basis(S),3) for S in sorted(subspaces)]
P=sum(1<<x for x in [0,1,2,3]);Q=sum(1<<x for x in [0,1,4,5]);T=sum(1<<x for x in [0,1,4,6,8,11])
large=[inspect([P,T],4),inspect([P,Q],4),inspect([T],4)]
assert [(x['seed_dimension'],x['closure_dimension']) for x in large]==[(14,15),(14,14),(15,16)]
assert sorted(map(len,large[1]['components']))==[1]*10+[6]
out={'scope':'all 67 seeds containing RM1 on 8 points, plus three exact 16-point models','small_seed_count':len(smallrows),'small_proper_closures':[r for r in smallrows if r['closure_dimension']>r['seed_dimension']],'large_models':large,'all_pair_partitions_equal':True,'python':platform.python_version(),'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
