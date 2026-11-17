#!/usr/bin/env python3
"""Exact two-seed completion DP; n=3 permutation cross-check and n=4 audit."""
import itertools,json
def run(n):
 N=n*n; full=(1<<N)-1; nb=[]
 for v in range(N):
  x,y=v%n,v//n; z=0
  for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
   X,Y=x+dx,y+dy
   if 0<=X<n and 0<=Y<n:z|=1<<(Y*n+X)
  nb.append(z)
 F=[0]*(1<<N);F[full]=1
 for mask in range(full-1,-1,-1):
  frontier=0
  q=mask
  while q:
   b=q&-q;v=b.bit_length()-1;frontier|=nb[v];q-=b
  frontier&=~mask
  while frontier:
   b=frontier&-frontier;F[mask]+=F[mask|b];frontier-=b
 def typ(v):
  x,y=v%n,v//n; return min((x,y),(n-1-x,y),(x,n-1-y),(n-1-x,n-1-y),(y,x),(n-1-y,x),(y,n-1-x),(n-1-y,n-1-x))
 rows=[]
 for p in range(N):
  for r in range(N):
   if p==r or nb[p]>>r&1:continue
   rows.append({'p':p,'r':r,'p_xy':[p%n,p//n],'r_xy':[r%n,r//n],'p_type':typ(p),'r_type':typ(r),'d':abs(p%n-r%n)+abs(p//n-r//n),'B':F[(1<<p)|(1<<r)]})
 return rows
r3=run(3);r4=run(4)
# same fixed p, and all moved roots in one D4 orbit
bad=[]
for p in range(16):
 for orbit in set(tuple(x['r_type']) for x in r4 if x['p']==p):
  a=sorted([x for x in r4 if x['p']==p and tuple(x['r_type'])==orbit],key=lambda x:x['d'])
  for x in a:
   for y in a:
    if x['d']<y['d'] and x['B']>y['B']:bad.append({'fixed_p':p,'root_orbit':orbit,'near':x,'far':y});break
   if bad:break
  if bad:break
 if bad:break
out={'n3_dp_sum':sum(x['B'] for x in r3),'n3_expected_B':81984,'n4_pairs':r4,'same_p_same_root_orbit_decreasing_example':bad[0] if bad else None,'scope':'exact subset DP; finite comparison only'}
assert out['n3_dp_sum']==81984
print(json.dumps(out,indent=2))
