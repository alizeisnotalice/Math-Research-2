#!/usr/bin/env python3
"""Independent 9! audit of usable early-interface directions."""
import itertools,json
from functools import lru_cache
from pathlib import Path
n=3;N=9;FULL=(1<<N)-1
adj=[0]*N
for y in range(n):
 for x in range(n):
  v=y*n+x
  for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
   X,Y=x+dx,y+dy
   if 0<=X<n and 0<=Y<n:adj[v]|=1<<(Y*n+X)
def comps(mask):
 out=[]
 while mask:
  seen=fr=mask&-mask
  while fr:
   z=fr&-fr;fr-=z;q=adj[z.bit_length()-1]&mask&~seen;seen|=q;fr|=q
  out.append(seen);mask&=~seen
 return out
@lru_cache(None)
def hdir(A,B,S=None):
 H=FULL^B
 if S is None:S=A
 if S==H:return 1
 return sum(hdir(A,B,S|(1<<v)) for v in range(N) if (H>>v&1) and not(S>>v&1) and (adj[v]&S))
cats={k:[0,0] for k in ('only_canonical_A','only_canonical_B','both')};two=0
for order in itertools.permutations(range(N)):
 seen=0;birth=[];merge=None
 for v in order:
  old=len(comps(seen)) if seen else 0
  if not(adj[v]&seen):birth.append(v)
  seen|=1<<v
  if merge is None and old==2 and len(comps(seen))==1:merge=(seen,v)
 if len(birth)!=2:continue
 two+=1;S,w=merge;A,B=comps(S^(1<<w));u0=hdir(A,B)>0;u1=hdir(B,A)>0
 assert u0 or u1
 key='both' if u0 and u1 else ('only_canonical_A' if u0 else 'only_canonical_B')
 d=(birth[0]%n-birth[1]%n)**2+(birth[0]//n-birth[1]//n)**2
 cats[key][0]+=1;cats[key][1]+=d
cpp=json.loads(Path('completion-n3.json').read_text())
assert two==143808
for k,(W,D) in cats.items():assert (W,D)==(int(cpp['availability'][k]['W']),int(cpp['availability'][k]['D']))
out={'status':'passed','n':3,'A2':two,'categories':{k:{'W':x[0],'D':x[1]} for k,x in cats.items()}}
Path('audit-completion-n3.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
