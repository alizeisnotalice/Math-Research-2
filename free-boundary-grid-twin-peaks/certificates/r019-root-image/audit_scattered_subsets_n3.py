#!/usr/bin/env python3
"""Exact 3x3 weighted injection for arbitrary shielded future subbranches."""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb
import json
from pathlib import Path

n=3;N=9
adj=[]
for y in range(n):
 for x in range(n):
  v=y*n+x;ns=[]
  for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
   X,Y=x+dx,y+dy
   if 0<=X<n and 0<=Y<n:ns.append(Y*n+X)
  adj.append(ns)
def peaks(word):
 pos=[0]*N
 for i,v in enumerate(word):pos[v]=i
 return [v for v in range(N) if all(pos[v]<pos[u] for u in adj[v])]
def connected(seq):
 if not seq:return False
 S=set(seq);seen={seq[0]};todo=[seq[0]]
 while todo:
  v=todo.pop()
  for u in adj[v]:
   if u in S and u not in seen:seen.add(u);todo.append(u)
 return len(seen)==len(S)
def growth(seq):return all(any(seq[i] in adj[seq[j]] for j in range(i)) for i in range(1,len(seq)))
def shuffles(a,b):
 for slots in combinations(range(len(a)+len(b)),len(a)):
  slots=set(slots);i=j=0;out=[]
  for k in range(len(a)+len(b)):
   if k in slots:out.append(a[i]);i+=1
   else:out.append(b[j]);j+=1
  yield tuple(out)

sources=events=invalid=0;weighted=Fraction();loads={};by_size={};zero=None
for word in permutations(range(N)):
 if len(peaks(word))!=1:continue
 sources+=1; local=Fraction()
 for a in range(1,N-1):
  A=word[:a];aset=set(A);w=word[a];future=word[a+1:];m=len(future)
  for sub in range(1,1<<m):
   B=tuple(future[i] for i in range(m) if sub>>i&1);b=len(B);bset=set(B)
   if min(A)>=min(B) or not connected(B) or not growth(B) or B[0] not in adj[w]:continue
   if any(y in adj[x] for x in aset for y in bset):continue
   tail=tuple(v for v in future if v not in bset); fibre=comb(m,b);gain=Fraction(comb(a+b,a),fibre)
   events+=1;weighted+=gain;local+=gain;by_size[b]=by_size.get(b,Fraction())+gain
   for prefix in shuffles(A,B):
    out=prefix+(w,)+tail
    if len(peaks(out))!=2:invalid+=1
    loads[out]=loads.get(out,Fraction())+Fraction(1,fibre)
 if local==0 and zero is None:zero=list(word)
maxload=max(loads.values(),default=Fraction())
if invalid or maxload>1:raise AssertionError((invalid,maxload))
out={'A1':sources,'events':events,'weighted_sum_num':weighted.numerator,'weighted_sum_den':weighted.denominator,
 'weighted_mean_num':weighted.numerator,'weighted_mean_den':weighted.denominator*sources,
 'weighted_image_mass_num':sum(loads.values(),Fraction()).numerator,'weighted_image_mass_den':sum(loads.values(),Fraction()).denominator,
 'distinct_outputs':len(loads),'max_load_num':maxload.numerator,'max_load_den':maxload.denominator,
 'invalid_outputs':invalid,'zero_example':zero,'by_B_size':{str(k):{'num':v.numerator,'den':v.denominator} for k,v in sorted(by_size.items())},'status':'passed'}
Path('scattered-subsets-n3.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
