"""Exact matching search for an MM completion in all 70 coordinate splits."""
import json,itertools,time,hashlib,sys
from pathlib import Path
P=Path(__file__).resolve().parent
INPUT=P.parent/'round032/exceptional_signed_mask_result.json'
r=json.loads(INPUT.read_text());S=r['points'];neg=r['solutions'][0]
physical={x:(neg>>i)&1 for i,x in enumerate(S)}
def fwht(a):
 a=a[:];h=1
 while h<len(a):
  for s in range(0,len(a),h*2):
   for i in range(s,s+h):a[i],a[i+h]=a[i]+a[i+h],a[i]-a[i+h]
  h*=2
 return a
A=[0]*256
for x,b in physical.items():A[x]=1-2*b
W=fwht(A);dual={a:int(w<0) for a,w in enumerate(W) if w}
assert len(dual)==32 and set(W)=={0,-16,16}
def proj(x,coords):return sum(((x>>c)&1)<<i for i,c in enumerate(coords))
def dot(x,y):return (x&y).bit_count()%2
def matching(adj):
 dest={}
 def dfs(y,seen):
  for a in adj[y]:
   if a in seen:continue
   seen.add(a)
   if a not in dest or dfs(dest[a],seen):dest[a]=y;return True
  return False
 for y in range(16):
  if not dfs(y,set()):return None
 return {y:a for a,y in dest.items()}
t=time.time();records=[];solution=None
for coords in itertools.combinations(range(8),4):
 other=tuple(i for i in range(8) if i not in coords)
 pc=[[] for _ in range(16)];dc=[[] for _ in range(16)]
 for x,b in physical.items():pc[proj(x,other)].append((proj(x,coords),b))
 for x,b in dual.items():dc[proj(x,coords)].append((proj(x,other),b))
 choices={};adj=[[] for _ in range(16)]
 for y in range(16):
  for a in range(16):
   for c in (0,1):
    if all(dot(x,a)^c==b for x,b in pc[y]) and all(dot(v,y)^c==b for v,b in dc[a]):
     choices[(y,a)]=c;adj[y].append(a);break
 pi=matching(adj);records.append(dict(data_coords=coords,edges=sum(map(len,adj)),perfect=pi is not None))
 if pi is not None:
  rho=[choices[y,pi[y]] for y in range(16)]
  f=[dot(proj(x,coords),pi[proj(x,other)])^rho[proj(x,other)] for x in range(256)]
  wh=fwht([1-2*b for b in f]);g=f[:]
  for x in S:g[x]^=1
  whg=fwht([1-2*b for b in g]);assert all(abs(w)==16 for w in wh+whg)
  solution=dict(data_coords=coords,selector_coords=other,pi=[pi[y] for y in range(16)],rho=rho,f_hex=hex(sum(b<<i for i,b in enumerate(f))),f_plus_h_hex=hex(sum(b<<i for i,b in enumerate(g))))
  break
out=dict(scope='70 coordinate MM splits, stop on first complete witness',splits_tested=len(records),records=records,solution=solution,seconds=time.time()-t,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),python=sys.version)
(P/'mm_completion_result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
