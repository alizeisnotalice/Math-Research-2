#!/usr/bin/env python3
"""Exact 3x3 audit of the low-peak label-gap lifting refinement."""
import itertools, json
from collections import Counter, defaultdict
from pathlib import Path
n=3; N=9
adj=[set() for _ in range(N)]
for y in range(n):
 for x in range(n):
  v=y*n+x
  for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
   X,Y=x+dx,y+dy
   if 0<=X<n and 0<=Y<n:adj[v].add(Y*n+X)
def dist(a,b):return abs(a%n-b%n)+abs(a//n-b//n)
def geo(w,d):
 return all(w[i+1] in adj[w[i]] for i in range(d-1)) and dist(w[0],w[d-1])==d-1
def fall(x,k):
 z=1
 for i in range(k):z*=x-i
 return z
a1=0; pref=Counter(); rows=defaultdict(Counter)
for w in itertools.permutations(range(N)):
 pos=[0]*N
 for i,v in enumerate(w):pos[v]=i
 peaks=[v for v in range(N) if all(pos[v]<pos[u] for u in adj[v])]
 if len(peaks)==1:
  a1+=1
  for d in range(2,5):
   if geo(w,d):pref[d]+=1
 elif len(peaks)==2:
  u,v=sorted(peaks,key=lambda x:pos[x]); d=dist(u,v)
  if 2<=d<=4:rows[d][N-1-pos[v]]+=1
out={'n':n,'A1':a1,'geodesic_prefix_counts':dict(pref),
     'distance_gap_counts':{str(d):dict(sorted(x.items())) for d,x in rows.items()},'checks':{}}
for d,h in rows.items():
 checks={}
 for g in range(N):
  lhs=sum(c for gap,c in h.items() if gap<=g)
  rhs=g*fall(N-1,d-2)*pref[d]
  assert lhs<=rhs
  checks[str(g)]={'A2_distance_d_low_peak_gap_le_g':lhs,'lifting_gap_upper_bound':rhs}
 out['checks'][str(d)]=checks
out['status']='passed'
Path('gap-audit-n3.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
