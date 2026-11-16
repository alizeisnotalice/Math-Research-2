#!/usr/bin/env python3
"""Independent 9! permutation audit of dominating-complement coverage."""
import itertools,json
from collections import Counter
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
  seen=front=mask&-mask
  while front:
   b=front&-front;front-=b;z=adj[b.bit_length()-1]&mask&~seen;seen|=z;front|=z
  out.append(seen);mask&=~seen
 return out
cats={k:[0,0] for k in ('only_A','only_B','both')};tw=0
for order in itertools.permutations(range(N)):
 seen=0;birth=[];merge=None
 for v in order:
  old=len(comps(seen)) if seen else 0
  if not(adj[v]&seen):birth.append(v)
  seen|=1<<v
  if merge is None and old==2 and len(comps(seen))==1:merge=(seen,v)
 if len(birth)!=2:continue
 tw+=1;S,w=merge;A,B=comps(S^(1<<w));C=FULL^(A|B)
 pa=all(adj[v]&A for v in range(N) if C>>v&1);pb=all(adj[v]&B for v in range(N) if C>>v&1)
 if not(pa or pb):continue
 key='both' if pa and pb else ('only_A' if pa else 'only_B')
 d=(birth[0]%n-birth[1]%n)**2+(birth[0]//n-birth[1]//n)**2
 cats[key][0]+=1;cats[key][1]+=d
cpp=json.loads(Path('dominating-n3.json').read_text())
assert tw==143808
for k,(w,d) in cats.items():assert (w,d)==(int(cpp['categories'][k]['W']),int(cpp['categories'][k]['D']))
out={'status':'passed','n':3,'A2':tw,'categories':{k:{'W':v[0],'D':v[1]} for k,v in cats.items()}}
Path('audit-dominating-n3.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
