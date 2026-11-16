#!/usr/bin/env python3
"""Exact 3x3 audit of the shielded scattered-lobe weighted injection."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import comb
import json
from pathlib import Path

n=3; N=9
adj=[]
for y in range(n):
 for x in range(n):
  v=y*n+x; q=[]
  for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
   X,Y=x+dx,y+dy
   if 0<=X<n and 0<=Y<n:q.append(Y*n+X)
  adj.append(tuple(q))
def peaks(word):
 pos=[0]*N
 for i,v in enumerate(word):pos[v]=i
 return [v for v in range(N) if all(pos[v]<pos[u] for u in adj[v])]
def components(vertices):
 unseen=set(vertices); ans=[]
 while unseen:
  seen={unseen.pop()}; todo=list(seen)
  while todo:
   v=todo.pop()
   for u in adj[v]:
    if u in unseen:unseen.remove(u);seen.add(u);todo.append(u)
  ans.append(seen)
 return ans
def shuffles(a,b):
 for slots in combinations(range(len(a)+len(b)),len(a)):
  slots=set(slots); i=j=0; out=[]
  for k in range(len(a)+len(b)):
   if k in slots:out.append(a[i]);i+=1
   else:out.append(b[j]);j+=1
  yield tuple(out)

sources=0; weighted=Fraction(); event_count=0; zero_K=None; loads={}; invalid=0
by_size=Counter()
for word in permutations(range(N)):
 if len(peaks(word))!=1:continue
 sources+=1; source_K=Fraction()
 for a in range(1,N-1):
  A=word[:a]; aset=set(A); w=word[a]; future=word[a+1:]
  for Bset in components(future):
   b=len(Bset)
   if min(A)>=min(Bset):continue
   if any(y in adj[x] for x in aset for y in Bset):continue
   beta=tuple(v for v in future if v in Bset)
   tail=tuple(v for v in future if v not in Bset)
   # Shielding forces the first beta vertex to attach to w; assert it explicitly.
   if beta[0] not in adj[w]:raise AssertionError(('bad root',word,a,Bset))
   c=N-a-b-1; fibre=comb(b+c,b); gain=Fraction(comb(a+b,a),fibre)
   source_K+=gain; weighted+=gain; event_count+=1;by_size[b]+=gain
   weight=Fraction(1,fibre)
   for prefix in shuffles(A,beta):
    out=prefix+(w,)+tail
    if len(peaks(out))!=2:invalid+=1
    loads[out]=loads.get(out,Fraction())+weight
 if source_K==0 and zero_K is None:zero_K=list(word)
maxload=max(loads.values(),default=Fraction())
if maxload>1:raise AssertionError(('fiber budget failed',maxload))
ans={
 'A1':sources,'events':event_count,
 'weighted_source_sum_numerator':weighted.numerator,'weighted_source_sum_denominator':weighted.denominator,
 'weighted_source_mean_numerator':weighted.numerator,'weighted_source_mean_denominator':weighted.denominator*sources,
 'weighted_image_mass_numerator':sum(loads.values(),Fraction()).numerator,
 'weighted_image_mass_denominator':sum(loads.values(),Fraction()).denominator,
 'distinct_outputs':len(loads),'max_weighted_output_load_numerator':maxload.numerator,'max_weighted_output_load_denominator':maxload.denominator,
 'invalid_outputs':invalid,'zero_weight_example':zero_K,
 'by_B_size':{str(k):{'numerator':v.numerator,'denominator':v.denominator} for k,v in sorted(by_size.items())},
 'status':'passed'
}
Path('scattered-lobes-n3.json').write_text(json.dumps(ans,indent=2)+'\n')
print(json.dumps(ans,indent=2))
