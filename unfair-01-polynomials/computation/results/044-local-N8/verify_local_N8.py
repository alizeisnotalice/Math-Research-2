#!/usr/bin/env python3
"""Rational interval certificate of a local-prefix counterexample, not original."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd,isqrt,lcm
import json,hashlib
D=Path(__file__).resolve().parent

def add(a,b):return [ (a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0)) for i in range(max(len(a),len(b))) ]
def scale(a,s):return [x*s for x in a]
def mul(a,b):
 r=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):r[i+j]+=x*y
 return r

def evalp(a,x):return sum(v*x**i for i,v in enumerate(a))
def rem(a,f):
 a=list(a)
 while len(a)>=len(f):
  c=a[-1]/f[-1];off=len(a)-len(f)
  for j,x in enumerate(f):a[off+j]-=c*x
  a.pop()
 return a

def main():
 q=F(59,64);aa=[[F(1)]]
 for k in range(1,9):
  a=add([F(1)],scale(aa[k-1],-q))
  if k>=4:a=add(a,[F(0)]+scale(aa[k-4],-1))
  aa.append(a)
 f=aa[8];assert len(f)==3 and f[-1]==1
 den=lcm(*(x.denominator for x in f));coeff=[int(x*den) for x in f];gg=gcd(*coeff);coeff=[x//gg for x in coeff]
 discr=coeff[1]**2-4*coeff[0]*coeff[2];assert discr>0 and isqrt(discr)**2!=discr
 lo,hi=F(4427,10000),F(4428,10000);assert evalp(f,lo)>0>evalp(f,hi);assert f[1]+2*hi<0
 intervals=[]
 for k,a in enumerate(aa[:8]):
  assert len(a)<=2
  vals=[evalp(a,lo),evalp(a,hi)];L,U=min(vals),max(vals)
  if k==0:assert L==U==1
  else:assert 0<L<=U<1
  intervals.append([str(L),str(U)])
 # B=(1,q,0,0,t); set A8=B8=0. Check full finite product in Q[t]/f.
 A=aa[:8]+[[F(0)]];B=[[F(1)],[q],[F(0)],[F(0)],[F(0),F(1)]]+[[F(0)]]*4
 C=[[F(0)] for _ in range(17)]
 for i,a in enumerate(A):
  for j,b in enumerate(B):C[i+j]=add(C[i+j],mul(a,b))
 for k in range(9):assert all(x==0 for x in rem(add(C[k],[-F(1)]),f))
 assert C[9]==[F(0)]+aa[5] # positive t*A5<1, so NOT binary beyond prefix.
 out=dict(status='passed_exact_local_N8_counterexample',q=str(q),t_minimal_polynomial_low_to_high=coeff,t_isolating_interval=[str(lo),str(hi)],A_polynomials_in_t=[[str(z) for z in a] for a in A],B_polynomials_in_t=[[str(z) for z in b] for b in B],A_intervals=intervals,verified_prefix=[0,8],A8=0,B8=0,first_fractional_pair=['5/64','59/64'],full_product_failure='C9=t*A5 lies strictly between 0 and 1',scope='Refutes only the all-N local prefix lemma; not the original normalized finite binary-product conjecture.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (D/'certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
