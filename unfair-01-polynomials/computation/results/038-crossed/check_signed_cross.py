"""Sparse integer verification; no symbolic factorizer or floating point."""
from pathlib import Path
from math import gcd
import json,os,hashlib,sys
assert os.path.ismount('/Volumes/data')
def mul(a,b):
 c={}
 for i,x in a.items():
  for j,y in b.items():c[i+j]=c.get(i+j,0)+x*y
 return {i:x for i,x in c.items() if x}
def sub(a,p):return {i*p:v for i,v in a.items()}
def binomials(weights):
 r={0:1}
 for w in weights:r=mul(r,{0:1,w:1})
 return r
def section(a,q):return {i//q:v for i,v in a.items() if i%q==0}
def binary(a):return set(a.values())=={1}
def family(p,q,b):
 assert 1<p<q and gcd(p,q)==1 and b>=3 and b%2
 F={p*q*k:(-1)**k for k in range(b)};G=binomials([p*p,p*q,q*q]);H=binomials([p*p,q*q,p*q*b])
 A=binomials([p*p*q*b,p*q*q,q**3]);B=binomials([p*q*q*b,p**3,p*p*q]);C=mul(A,B)
 assert mul(F,G)==H and mul(sub(F,p),sub(G,q))==A and mul(sub(F,q),sub(G,p))==B
 assert C==mul(sub(H,p),sub(H,q));assert all(binary(v) for v in [G,H,A,B,C])
 assert F[p*q]==-1 and F[0]==F[max(F)]==1
 return F,G,H,A,B,C
cases=[]
for q in range(3,18):
 for p in range(2,q):
  if gcd(p,q)>1:continue
  S=(p+q)*(p*p+q*q);b=3
  while p*q*b*min(p,q-p)<=S:b+=2
  F,G,H,A,B,C=family(p,q,b)
  RF,RG,H0=section(F,q),section(G,q),section(H,q)
  U=mul(F,sub(RG,p));V=mul(G,sub(RF,p))
  assert mul(RF,RG)==H0 and U==section(B,q) and V==section(A,q)
  assert binary(mul(U,V)) and max(mul(U,V))<max(C)
  cases.append({'p':p,'q':q,'b':b,'degrees':[max(F),max(G),max(H),max(C)],'C_terms':len(C)})
# Compact full-model counterexample to deleting the nonzero residue carry.
p,q,b,L,M=2,3,19,16,446
F,G,H,A,B,C=family(p,q,b)
F=mul(F,binomials([L]));G=mul(G,binomials([M]));H=mul(H,binomials([L,M]))
A=mul(A,binomials([p*L,q*M]));B=mul(B,binomials([q*L,p*M]));C=mul(A,B)
assert mul(F,G)==H and mul(sub(F,p),sub(G,q))==A and mul(sub(F,q),sub(G,p))==B
assert C==mul(sub(H,p),sub(H,q));assert all(binary(v) for v in [G,H,A,B,C]);assert F[6]==-1
assert min(max(F),max(G))>=2*q
RF,RG,H0=section(F,q),section(G,q),section(H,q);S=mul(RF,RG)
assert RF=={2*k:(-1)**k for k in range(19)}
assert RG==binomials([2,3,150]) and H0==mul(binomials([38,3]),{0:1,150:1,154:1})
assert H0!=S
carry={k:H0.get(k,0)-S.get(k,0) for k in set(H0)|set(S)};carry={k:v for k,v in carry.items() if v}
assert carry=={154+k:v for k,v in binomials([38,3]).items()}
U,V=mul(F,sub(RG,p)),mul(G,sub(RF,p));UV=mul(U,V);C0=section(C,q)
assert U==section(B,q) and V==section(A,q) and binary(UV)
assert all(0<=UV.get(k,0)<=C0.get(k,0) for k in set(UV)|set(C0))
record={'status':'passed','scope':'auxiliary maternal-positivity and zero-carry counterexamples, not unfair original factors','python':sys.version,'cases':cases,'odd_example':next(c for c in cases if(c['p'],c['q'])==(3,5)), 'carry_example':{'p':p,'q':q,'b':b,'L':L,'M':M,'degrees':[max(F),max(G),max(H),max(C)],'support_sizes':[len(F),len(G),len(H),len(A),len(B),len(C)],'carry':carry,'projected_product_terms':len(UV),'C_zero_section_terms':len(C0),'projected_product_degree':max(UV)},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k!='cases'}));print('finite sanity cases',len(cases))
