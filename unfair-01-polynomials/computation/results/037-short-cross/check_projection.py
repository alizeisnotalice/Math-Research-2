from fractions import Fraction as Q
import itertools,hashlib,json,sys,os
from pathlib import Path
assert os.path.ismount('/Volumes/data')
def mul(a,b):
 c=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def sub(a,p):
 c=[Q(0)]*((len(a)-1)*p+1)
 for i,v in enumerate(a):c[i*p]=v
 return c
def trim(a):
 while len(a)>1 and not a[-1]:a.pop()
 return a
def project(a,q):return trim(a[::q])
count=0
for p,q in [(2,3),(2,21),(3,22),(5,23)]:
 for coeffs in itertools.product([Q(-1),Q(0),Q(1,2),Q(1)],repeat=3):
  F=[Q(1),*coeffs,Q(1)];G=[Q(0)]*(2*q)
  G[0]=G[q]=G[-1]=Q(1);G[q-1]=Q(2,3)
  B=mul(sub(F,q),sub(G,p));R=[Q(1)]+[Q(0)]*(p-1)+[Q(1)]
  assert project(B,q)==trim(mul(F,R));count+=1
# Legal full model: (1+x^2)(1+x^21) has distinct sums and crossed product binary.
p,q=2,21;F=[Q(1),Q(1)];G=[Q(1)]+[Q(0)]*2+[Q(1)]
H=mul(F,G);C=mul(sub(H,p),sub(H,q));assert all(x in(0,1) for x in C)
# Critical invalid generalization: adding g_(2q) adds a third coefficient to the extracted slice.
F=[Q(1),Q(-1),Q(1)];G=[Q(0)]*(2*q+1);G[0]=G[q]=G[2*q]=Q(1)
actual=project(mul(sub(F,q),sub(G,p)),q)
R=[Q(1)]+[Q(0)]*(p-1)+[Q(1)]
assert actual!=trim(mul(F,R))
# No-gap gate: H=1+x^p violates it and (1+x^p)H has coefficient2.
assert max(mul(R,R))==2
out={'status':'passed','exact_projection_cases':count,'negative_regressions':['degree_equal_2q_adds_term','missing_no_gap_creates_coefficient2'],'python':sys.version,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'coefficient identity and explicit gate checks; analytic proof supplies positivity and degree descent'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
