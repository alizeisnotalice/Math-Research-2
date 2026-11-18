#!/usr/bin/env python3
"""Exact necessary rows 5--8 continuation for degree-8 S={2}.

Input is the frozen stage-4 ledger.  This script keeps G4=Q4=0 and forces
A6=B6=0, but deliberately does *not* set G6 or Q6 to zero.
"""
from __future__ import annotations
import argparse, concurrent.futures as cf, hashlib, json, math, os, sys, time
from collections import Counter
from pathlib import Path

Pair=tuple[int,int]; ZERO=(0,0); ONE=(1,0)
BOUNDS={5:(1025,8),6:(4097,16),7:(16385,32),8:(65537,64)}

def add(x,y):return (x[0]+y[0],x[1]+y[1])
def scale(k,x):return (k*x[0],k*x[1])
def mul(x,y,T,N):
 p,q=x;r,s=y;return (p*r-q*s*N,p*s+q*r+q*s*T)
def neg(x):return (-x[0],-x[1])
def conj(x,T):p,q=x;return (p+q*T,-q)
def sign_plus(x,T,D):
 p,q=x;r=2*p+q*T
 if q==0:return (r>0)-(r<0)
 if q>0 and r>=0:return 1
 if q<0 and r<=0:return -1
 z=(r*r>q*q*D)-(r*r<q*q*D)
 return -z if q>0 else z
def floor_plus(x,T,D):
 p,q=x;R=abs(p)+abs(q)*(abs(T)+math.isqrt(D)+3)+4;lo,hi=-R,R
 while hi-lo>1:
  m=(lo+hi)//2
  if sign_plus((p-m,q),T,D)>=0:lo=m
  else:hi=m
 return lo
def shifts(base,T,D):
 if base[1]==0:return [(-base[0],ZERO),(1-base[0],ONE)]
 h=-floor_plus(base,T,D);z=(base[0]+h,base[1])
 assert sign_plus(z,T,D)>0 and sign_plus(add(ONE,neg(z)),T,D)>0
 return [(h,z)]
def binary(x):return x[1]==0 and x[0] in (0,1)
def both_abs(x,bound,T,D):
 num,den=bound;r=abs(2*x[0]+x[1]*T);q=abs(x[1]);gap=2*num-r*den
 return gap>0 and q*q*D*den*den<gap*gap
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.partial')
 with t.open('w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)

def newton(L,T,N,b,c,d,e,f):
 co=[None,(L,0),(0,1),(b,0),(c,0),(d,0),(e,0),(f,0),(1,0)]
 s={0:(8,0)}
 for k in range(1,9):
  z=ZERO
  for j in range(1,k):z=add(z,mul(co[j],s[k-j],T,N))
  z=add(z,scale(k,co[k]));s[k]=neg(z)
 return s

def symbolic_check():
 import sympy as sp
 L,a,b,c,d,e,f=sp.symbols('L a b c d e f'); G1,G2,G3,G4,G5,G6,G7=sp.symbols('G1:8')
 P=[1,L,a,b,c,d,e,f,1];G=[1,G1,G2,G3,G4,G5,G6,G7]
 A=[sp.expand(sum(P[i]*G[k-i] for i in range(k+1))) for k in range(8)]
 assert sp.expand(A[5]-(G5+d+a*G3+b*G2+c*G1+L*G4))==0
 assert sp.expand(A[6]-(G6+e+L*G5+a*G4+b*G3+c*G2+d*G1))==0
 assert sp.expand(A[7]-(G7+f+L*G6+a*G5+b*G4+c*G3+d*G2+e*G1))==0
 assert sp.diff(A[6],a)==G4
 T,N=sp.symbols('T N'); ss={0:8}
 co=[None,L,a,b,c,d,e,f,1]
 for k in range(1,9):ss[k]=sp.expand(-sum(co[j]*ss[k-j] for j in range(1,k))-k*co[k])
 # Substitute a^2 recursively just to independently check the recurrence
 # is present through p8; exact pair arithmetic below performs reduction.
 assert ss[5].has(d) and ss[6].has(e) and ss[7].has(f) and ss[8].has(sp.Integer(-8))

def load_stage4(path:Path):
 s=json.loads((path/'summary.json').read_text())
 assert s['status']=='complete' and s['cores']==s['completed']==11760 
 fs=sorted((path/'cores').glob('*.json'));assert len(fs)==11760
 rawhash=hashlib.sha256();jobs=[];total=0
 for p in fs:
  raw=p.read_bytes();rawhash.update(p.name.encode()+b'\0'+raw);x=json.loads(raw)
  assert x['status']=='complete';total+=len(x['records']); jobs.append((p,sha(p)))
 assert total==s['row4_states']
 return jobs,hashlib.sha256((sha(path/'summary.json')+rawhash.hexdigest()).encode()).hexdigest()

def one(job):
 source,source_hash,outdir,input_hash,codehash,maxseconds=job; x=json.loads(Path(source).read_text());L,T,N,D=x['binding']['core'];begin=time.monotonic();ct=Counter();final=[]
 target=Path(outdir)/'cores'/Path(source).name
 bind={'version':1,'core':[L,T,N,D],'stage4_source_sha256':source_hash,'stage4_input_sha256':input_hash,'code_sha256':codehash,'bounds':{'d':[-92,92],'e':[-41,41],'f':[-9,9],'p5':[1025,8],'p6':[4097,16],'p7':[16385,32],'p8':[65537,64]}}
 if target.exists():
  o=json.loads(target.read_text())
  if o.get('status')=='complete' and o.get('binding')==bind:return o['summary']
  raise RuntimeError(f'stale cache {target}')
 for r in x['records']:
  rL,rT,rN,rD,b,c,G1,Q1,rho,sigma,eta,zeta=r
  assert (rL,rT,rN,rD)==(L,T,N,D)
  A1=(L+G1,0);B1=(L+Q1,0);A2=(rho,1);B2=conj((sigma,1),T)
  G2=rho-L*G1;Q2=sigma-L*Q1;G3=eta-b;Q3=zeta-b
  U3=(L*G2,G1);V3=add((L*Q2,0),conj((0,Q1),T));A3=add(U3,(eta,0));B3=add(V3,(zeta,0))
  W4=(L*G3+b*G1,G2);Z4=add((L*Q3+b*Q1,0),conj((0,Q2),T));A4=add(W4,(c,0));B4=add(Z4,(c,0))
  assert all((binary(z) or z[1]!=0) for z in (A1,B1,A2,B2,A3,B3,A4,B4))
  # G4=Q4=0; U5 is independent of d, and row-five unit choices retain
  # both endpoints when a base happens to be rational.
  U5=(b*G2+c*G1,G3);V5=add((b*Q2+c*Q1,0),conj((0,Q3),T))
  for eta5,A5 in shifts(U5,T,D):
   for zeta5,B5 in shifts(V5,T,D):
    C5=add(add(A5,B5),add(add(scale(A1[0],B4),mul(A2,B3,T,N)),add(mul(A3,B2,T,N),scale(B1[0],A4))))
    if not binary(C5):ct['c5_reject']+=1;continue
    ct['row5_states']+=1
    C6=ZERO
    Avec=[ONE,A1,A2,A3,A4,A5,ZERO];Bvec=[ONE,B1,B2,B3,B4,B5,ZERO]
    for i in range(1,6):C6=add(C6,mul(Avec[i],Bvec[6-i],T,N))
    if not binary(C6):ct['c6_reject']+=1;continue
    ct['cross6_pass']+=1
    for d in range(-92,93):
     ct['d_box']+=1;s=newton(L,T,N,b,c,d,0,0)
     if not both_abs(s[5],BOUNDS[5],T,D):ct['s5_reject']+=1;continue
     ct['s5_pass']+=1;G5=eta5-d;Q5=zeta5-d
     for e in range(-41,42):
      ct['e_box']+=1;s=newton(L,T,N,b,c,d,e,0)
      if not both_abs(s[6],BOUNDS[6],T,D):ct['s6_reject']+=1;continue
      ct['s6_pass']+=1
      G6=-L*G5-b*G3-c*G2-d*G1-e;Q6=-L*Q5-b*Q3-c*Q2-d*Q1-e
      assert add((L*G5+b*G3+c*G2+d*G1+e+G6,0),ZERO)==ZERO
      U7=(L*G6+c*G3+d*G2+e*G1,G5);V7=add((L*Q6+c*Q3+d*Q2+e*Q1,0),conj((0,Q5),T))
      for eta7,A7 in shifts(U7,T,D):
       for zeta7,B7 in shifts(V7,T,D):
        C7=add(A7,B7)
        Avec=[ONE,A1,A2,A3,A4,A5,ZERO,A7];Bvec=[ONE,B1,B2,B3,B4,B5,ZERO,B7]
        for i in range(1,7):C7=add(C7,mul(Avec[i],Bvec[7-i],T,N))
        if not binary(C7):ct['c7_reject']+=1;continue
        ct['row7_states']+=1
        for f in range(-9,10):
         ct['f_box']+=1;s=newton(L,T,N,b,c,d,e,f)
         if not both_abs(s[7],BOUNDS[7],T,D):ct['s7_reject']+=1;continue
         if not both_abs(s[8],BOUNDS[8],T,D):ct['s8_reject']+=1;continue
         ct['row8_states']+=1
         final.append([L,T,N,D,b,c,d,e,f,G1,Q1,rho,sigma,eta,zeta,eta5,zeta5,G6,Q6,eta7,zeta7])
  if time.monotonic()-begin>maxseconds:raise TimeoutError(f'core cap {L,T,N,D}')
 summary={'L':L,'T':T,'N':N,'Delta':D,'counts':dict(ct),'row8_states':len(final),'seconds':time.monotonic()-begin}
 atomic(target,{'status':'complete','binding':bind,'summary':summary,'records':final});return summary

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stage4',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=1);ap.add_argument('--max-seconds',type=float,default=14400);ap.add_argument('--limit',type=int);a=ap.parse_args()
 symbolic_check();jobs,input_hash=load_stage4(a.stage4)
 if a.limit is not None:jobs=jobs[:a.limit]
 codehash=sha(Path(__file__));work=[(*q,str(a.out),input_hash,codehash,a.max_seconds) for q in jobs];res=[]
 with cf.ProcessPoolExecutor(max_workers=a.workers) as ex:
  for z in ex.map(one,work):res.append(z)
 ct=Counter();[ct.update(z['counts']) for z in res]
 atomic(a.out/'summary.json',{'status':'complete','scope':'Necessary degree-8 S={2} rows5--8 from frozen stage4 only; no root screen or arbitrary-cofactor conclusion.','code_sha256':codehash,'stage4_input_sha256':input_hash,'cores':len(jobs),'completed':len(res),'counts':dict(ct),'row8_states':sum(z['row8_states'] for z in res),'workers':a.workers,'python':sys.version})
 print(json.dumps({'status':'complete','cores':len(jobs),'row8_states':sum(z['row8_states'] for z in res),'counts':dict(ct)},sort_keys=True))
if __name__=='__main__':main()
