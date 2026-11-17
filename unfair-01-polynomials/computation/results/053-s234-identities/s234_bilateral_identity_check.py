#!/usr/bin/env python3
"""Independent bilateral formal check of the S={2,3,4} third-row relations."""
import hashlib,json,os
from collections import Counter
from pathlib import Path
OUT=Path(__file__).with_suffix('.json')
V=('L','A','Ah','K','g2','q2'); Z=(0,)*len(V)
def atomvar(name):
 e=list(Z);e[V.index(name)]=1;return {tuple(e):1}
def one():return {Z:1}
def add(*xs):
 o=Counter()
 for x in xs:o.update(x)
 return {e:c for e,c in o.items() if c}
def neg(x):return {e:-c for e,c in x.items()}
def mul(x,y):
 o=Counter()
 for e,c in x.items():
  for f,d in y.items():o[tuple(u+v for u,v in zip(e,f))]+=c*d
 return {e:c for e,c in o.items() if c}
L,A,Ah,K,g2,q2=map(atomvar,V)
P=[one(),L,A,add(mul(L,A),K)]
Pb=[one(),L,Ah,add(mul(L,Ah),K)]
G=[one(),neg(L),g2,neg(add(mul(L,g2),K))]
Q=[one(),neg(L),q2,neg(add(mul(L,q2),K))]
def cv(X,Y,n):return add(*(mul(X[i],Y[n-i]) for i in range(n+1)))
A1,A2,A3=[cv(G,P,n) for n in (1,2,3)]
B1,B2,B3=[cv(Q,Pb,n) for n in (1,2,3)]
assert A1=={} and B1=={} and A3=={} and B3=={}
assert A2==add(A,g2,neg(mul(L,L)))
assert B2==add(Ah,q2,neg(mul(L,L)))
obj={'status':'passed','checks':['G1=Q1=-ell','G3=-ell*G2-k','Q3=-ell*Q2-k','A1=B1=A3=B3=0','A2=a+G2-ell^2','B2=abar+Q2-ell^2'],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
t=OUT.with_suffix('.partial')
with t.open('w') as f:json.dump(obj,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
os.replace(t,OUT);print(json.dumps(obj,sort_keys=True))
