#!/usr/bin/env python3
"""Independent exact verifier of the degree-8 S={1} root-screen ledger.

It does not import the production p+q*a code.  It uses
P=F+a*x and H=F^2+T*x*F+N*x^2 over Z[x], then repeats positive-axis and
two scaled Schur-Cohn necessary tests for every deduplicated mother.
"""
import argparse,hashlib,json,math,os
from collections import Counter
from pathlib import Path

POINTS=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def conv(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out
def H(row):
 T,N,b,c,d,e,f,g=row;F=[1,0,b,c,d,e,f,g,1];h=conv(F,F)
 for i,x in enumerate(F):h[i+1]+=T*x
 h[2]+=N
 assert len(h)==17 and h[0]==h[-1]==1
 return h
def schur(co):
 """Strict Schur-Cohn descent, coefficients from highest to constant."""
 a=list(co)
 while len(a)>1:
  d=math.gcd(*a)
  if a[0]<0:d=-d
  a=[x//d for x in a]
  if abs(a[-1])>=a[0]:return len(a)-1
  lead,constant=a[0],a[-1];degree=len(a)-1
  a=[lead*a[i]-constant*a[degree-i] for i in range(degree)]
 return None
def screen(row):
 h=H(row)
 for j,(p,q) in enumerate(POINTS):
  if sum(h[i]*p**i*q**(16-i) for i in range(17))<=0:return 10+j
 k=schur([h[i]*(2**i) for i in range(16,-1,-1)])
 if k is not None:return 30+k
 k=schur([h[i]*(2**(16-i)) for i in range(17)])
 if k is not None:return 50+k
 return 255
def atomic(p,x):
 t=p.with_suffix(p.suffix+'.partial')
 with t.open('w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stage57',type=Path,required=True);ap.add_argument('--production',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 rows=set()
 for p in sorted((a.stage57/'cores').glob('*.json')):
  x=json.loads(p.read_text())
  for b,c,d,e,f,g,*_ in x['records']:
   q=x['summary'];rows.add((q['T'],q['N'],b,c,d,e,f,g))
 assert len(rows)==67
 expected={tuple(x['binding']['row']):x['screen']['code'] for p in a.production.glob('mother-*.json') for x in [json.loads(p.read_text())]}
 assert set(expected)==rows and len(expected)==67
 got={r:screen(r) for r in sorted(rows)}
 assert got==expected
 out={'status':'passed_all_67','stage57_sha256':sha(a.stage57/'summary.json'),'production_summary_sha256':sha(a.production/'summary.json'),'rows':67,'codes':dict(Counter(got.values())),'source_sha256':sha(Path(__file__))}
 atomic(a.out,out);print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
