#!/usr/bin/env python3
"""Exact continuum prefix verification by rational Bernstein coefficients.
No discovery script or symbolic algebra package is imported.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import argparse,hashlib,json,os,sys,time

def trim(v):
 v=list(v)
 while len(v)>1 and not v[-1]:v.pop()
 return v

def add(a,b):
 c=[Q(0)]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return trim(c)

def neg(a):return [-x for x in a]
def mul(a,b):
 c=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)
def conv(A,B):
 out={}
 for i,a in A.items():
  for j,b in B.items():out[i+j]=add(out.get(i+j,[Q(0)]),mul(a,b))
 return out

def bernstein(poly,L):
 n=len(poly)-1
 return [sum((poly[i]*L**i*Q(comb(k,i),comb(n,i)) for i in range(k+1)),Q(0)) for k in range(n+1)]

def verify(source,s,L):
 H=source['first']['H_prefix'];assert len(H)==42 and all(type(v)is int and v in (0,1) for v in H)
 G={0:[Q(1)],s:[Q(0),Q(1)],20:[Q(1)],21:[Q(1)],21+s:[Q(0),Q(1)],41:[Q(1)]}
 F={}
 for i in range(42):
  v=[Q(H[i])]
  for j,g in G.items():
   if 0<j<=i:v=add(v,neg(mul(g,F[i-j])))
  F[i]=v
 product=conv(F,G)
 assert all(product.get(i,[Q(0)])==[Q(H[i])] for i in range(42))
 constraints=[('F',i,v) for i,v in F.items()]
 for name,u,v,limit in [('A',2,21,84),('B',21,2,882)]:
  C=conv({i*u:x for i,x in F.items()},{j*v:y for j,y in G.items()})
  for i in range(limit):
   z=C.get(i,[Q(0)]);constraints.extend([(name+'lower',i,z),(name+'upper',i,add([Q(1)],neg(z)))])
 failures=[];unique={}
 for name,i,poly in constraints:
  key=tuple(poly);bs=bernstein(poly,L)
  if any(x<0 for x in bs):failures.append((name,i,list(map(str,poly)),list(map(str,bs))))
  unique[key]=bs
 assert not failures,failures
 # Literal binary C on the entire chosen H-prefix polynomial.
 hh={i:[Q(v)] for i,v in enumerate(H) if v}
 C=conv({2*i:v for i,v in hh.items()},{21*i:v for i,v in hh.items()})
 assert all(v==[Q(1)] for v in C.values())
 if s==8:assert F[36]==[Q(0),Q(1),Q(-2)]
 rational=[sum((c*Q(1,2)**j for j,c in enumerate(F[i])),Q(0)) for i in range(42)]
 assert rational==list(map(Q,source['first']['F_prefix']))
 return {'s':s,'closed_interval':[0,str(L)],'constraints':len(constraints),'distinct_polynomials':len(unique),'F_nonzero_polynomials':{str(i):list(map(str,v)) for i,v in F.items() if v!=[Q(0)]},'bernstein_bounds':[{'power_coefficients':list(map(str,key)),'bernstein_coefficients':list(map(str,val))} for key,val in sorted(unique.items())],'necessity_for_s8':'F36=a(1-2a); under 0<a<1, F36>=0 iff a<=1/2' if s==8 else 'domain 0<a<1 assumed; all of it certified'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--forward',type=Path,required=True);ap.add_argument('--reverse',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert os.path.ismount('/Volumes/data');t=time.monotonic()
 data=[json.loads(a.forward.read_text()),json.loads(a.reverse.read_text())]
 result={'scope':'separate finite forward/reciprocal prefix witnesses; no single complete F is constructed','results':[verify(data[0],8,Q(1,2)),verify(data[1],12,Q(1))],'python':sys.version,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':[hashlib.sha256(p.read_bytes()).hexdigest() for p in [a.forward,a.reverse]],'elapsed_seconds':time.monotonic()-t}
 a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','elapsed_seconds':result['elapsed_seconds'],'distinct_constraints':[r['distinct_polynomials'] for r in result['results']]}))
if __name__=='__main__':main()
