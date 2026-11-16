#!/usr/bin/env python3
from itertools import product
from fractions import Fraction as F
from bisect import bisect_right
from pathlib import Path
import json,sys
checked=0;uniform=0;strict=0
for P in range(1,7):
 for a in product(range(4),repeat=P):
  T=sum(a)
  if T==0:continue
  c=[];s=0
  for v in a:s+=v;c.append(s)
  counts=[0]*P;var=F();means=[]
  for j in range(P):
   g=[]
   for t in range(T):
    u=(j*T+t)//P;i=bisect_right(c,u);assert 0<=u<T and a[i]>0
    counts[i]+=1;g.append(i*i)
    if len(set(a))==1:assert i==j
   mu=F(sum(g),T);means.append(mu);var+=(F(sum(x*x for x in g),T)-mu*mu)/(P*P)
  assert counts==[P*v for v in a]
  mu=sum(F(a[i]*i*i,T) for i in range(P));v=sum(F(a[i]*i**4,T) for i in range(P))-mu*mu
  assert sum(means)/P==mu and 0<=var<=v/P
  strict+=var<v/P;uniform+=len(set(a))==1;checked+=1
out={'status':'passed','weight_vectors':checked,'uniform_weight_vectors':uniform,'strict_variance_improvements_for_test_i_squared':strict,'scope':'finite exact verification; general lemma is proved separately'}
p=Path(sys.argv[1]) if len(sys.argv)>1 else None
if p:p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
