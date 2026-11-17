#!/usr/bin/env python3
"""Exact necessary rows 5--7 from frozen S={1}, degree-eight stage4 prefixes.

Independent of the stage4 implementation: it reads its JSON records but
rebuilds all raw G/Q rows and all A/B/C coefficients itself.  No finite-depth
survival is interpreted as a candidate or an exclusion.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path

Pair=tuple[int,int]  # p+q*a, where a^2=T*a-N
ZERO=(0,0); ONE=(1,0)
# Lists (rather than tuples) make checkpoint bindings stable across JSON round trips.
BOUNDS={'5':[1025,8],'6':[4097,16],'7':[16385,32],'8':[65537,64]}

def add(x:Pair,y:Pair)->Pair:return (x[0]+y[0],x[1]+y[1])
def neg(x:Pair)->Pair:return (-x[0],-x[1])
def scale(k:int,x:Pair)->Pair:return (k*x[0],k*x[1])
def mul(x:Pair,y:Pair,T:int,N:int)->Pair:
 p,q=x;r,s=y;return (p*r-q*s*N,p*s+q*r+q*s*T)
def conj(x:Pair,T:int)->Pair:
 p,q=x;return (p+q*T,-q)

def sign(x:Pair,T:int,D:int)->int:
 p,q=x;r=2*p+q*T
 if q==0:return (r>0)-(r<0)
 if r>=0 and q>0:return 1
 if r<=0 and q<0:return -1
 left=r*r;right=q*q*D
 if r>0:return (left>right)-(left<right)
 return (right>left)-(right<left)

def floor_pair(x:Pair,T:int,D:int)->int:
 p,q=x;B=abs(p)+abs(q)*(abs(T)+math.isqrt(D)+2)+3;lo,hi=-B,B
 while hi-lo>1:
  m=(lo+hi)//2
  if sign((p-m,q),T,D)>=0:lo=m
  else:hi=m
 return lo

def shifts(base:Pair,T:int,D:int)->list[tuple[int,Pair]]:
 p,q=base
 if q==0:return [(-p,ZERO),(1-p,ONE)]
 h=-floor_pair(base,T,D);z=(p+h,q)
 assert sign(z,T,D)>0 and sign(add(ONE,neg(z)),T,D)>0
 return [(h,z)]

def within(x:Pair,num:int,den:int,T:int,D:int)->bool:
 p,q=x
 return sign((num-den*p,-den*q),T,D)>0 and sign((num+den*p,den*q),T,D)>0
def both(x:Pair,num:int,den:int,T:int,D:int)->bool:
 return within(x,num,den,T,D) and within(conj(x,T),num,den,T,D)

def binary(A:tuple[Pair,...],B:tuple[Pair,...],k:int,T:int,N:int)->int|None:
 z=ZERO
 for i in range(k+1):z=add(z,mul(A[i],B[k-i],T,N))
 return z[0] if z[1]==0 and z[0] in (0,1) else None

def newton_base4(a:Pair,b:int,c:int,d:int,T:int,N:int)->tuple[Pair,Pair,Pair,Pair]:
 s1=neg(a)
 s2=add(neg(mul(a,s1,T,N)),(-2*b,0))
 s3=add(add(neg(mul(a,s2,T,N)),neg(scale(b,s1))),(-3*c,0))
 s4=add(add(add(neg(mul(a,s3,T,N)),neg(scale(b,s2))),neg(scale(c,s1))),(-4*d,0))
 return s1,s2,s3,s4

def atomic(path:Path,obj:object)->None:
 t=path.with_suffix(path.suffix+'.partial')
 with t.open('w') as f:json.dump(obj,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,path)
def sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()

def one(src:str,srcsha:str,outdir:str,codehash:str,maxseconds:float)->dict:
 source=Path(src);raw=source.read_bytes();assert hashlib.sha256(raw).hexdigest()==srcsha
 old=json.loads(raw);assert old['status']=='complete'
 T,N,D=old['summary']['T'],old['summary']['N'],old['summary']['Delta']
 binding={'code_sha256':codehash,'stage4_file_sha256':srcsha,'core':[T,N,D],'bounds':BOUNDS}
 target=Path(outdir)/'cores'/source.name
 if target.exists():
  z=json.loads(target.read_text());assert z['binding']==binding and z['status']=='complete';return z['summary']|{'cached':True}
 begin=time.monotonic();a=(0,1);abar=conj(a,T);g1=-floor_pair(a,T,D);q1=-floor_pair(abar,T,D)
 A1=(g1,1);B1=(q1+T,-1)
 assert binary((ONE,A1),(ONE,B1),1,T,N)==1
 ct=Counter(stage4_records=len(old['records']));final=[]
 for b,c,d,rho,sigma,eta,zeta,alpha,beta in old['records']:
  G2=rho-b;Q2=sigma-b;G3=eta-c;Q3=zeta-c;G4=alpha-d-b*rho+b*b;Q4=beta-d-b*sigma+b*b
  A2=add(scale(g1,a),(rho,0));B2=add(scale(q1,abar),(sigma,0))
  A3=add(add(scale(rho,a),scale(b,add((g1,0),neg(a)))),(eta,0));B3=add(add(scale(sigma,abar),scale(b,add((q1,0),neg(abar)))),(zeta,0))
  A4=add(add(scale(eta,a),scale(c,add((g1,0),neg(a)))),(alpha,0));B4=add(add(scale(zeta,abar),scale(c,add((q1,0),neg(abar)))),(beta,0))
  assert binary((ONE,A1,A2),(ONE,B1,B2),2,T,N) is not None
  assert binary((ONE,A1,A2,A3),(ONE,B1,B2,B3),3,T,N) is not None
  assert binary((ONE,A1,A2,A3,A4),(ONE,B1,B2,B3,B4),4,T,N) is not None
  s1,s2,s3,s4=newton_base4(a,b,c,d,T,N)
  U5=add(add(add(scale(G4,a),scale(b,(G3,0))),scale(c,(G2,0))),scale(d,(g1,0)))
  V5=add(add(add(scale(Q4,abar),scale(b,(Q3,0))),scale(c,(Q2,0))),scale(d,(q1,0)))
  for eta5,A5 in shifts(U5,T,D):
   for zeta5,B5 in shifts(V5,T,D):
    if binary((ONE,A1,A2,A3,A4,A5),(ONE,B1,B2,B3,B4,B5),5,T,N) is None:ct['c5_reject']+=1;continue
    ct['row5_states']+=1
    for e in range(-92,93):
     ct['e_box']+=1
     s5=add(add(add(add(neg(mul(a,s4,T,N)),neg(scale(b,s3))),neg(scale(c,s2))),neg(scale(d,s1))),(-5*e,0))
     if not both(s5,*BOUNDS['5'],T,D):ct['s5_reject']+=1;continue
     ct['s5_pass']+=1
     G5=eta5-e;Q5=zeta5-e
     U6=add(add(add(add(scale(G5,a),scale(b,(G4,0))),scale(c,(G3,0))),scale(d,(G2,0))),scale(e,(g1,0)))
     V6=add(add(add(add(scale(Q5,abar),scale(b,(Q4,0))),scale(c,(Q3,0))),scale(d,(Q2,0))),scale(e,(q1,0)))
     su={h for h,_ in shifts(U6,T,D)};sv={h for h,_ in shifts(V6,T,D)}
     for f in sorted(su&sv&set(range(-41,42))):
      A6=add(U6,(f,0));B6=add(V6,(f,0));assert sign(A6,T,D)>=0 and sign(add(ONE,neg(A6)),T,D)>=0
      assert sign(B6,T,D)>=0 and sign(add(ONE,neg(B6)),T,D)>=0
      s6=add(add(add(add(add(neg(mul(a,s5,T,N)),neg(scale(b,s4))),neg(scale(c,s3))),neg(scale(d,s2))),neg(scale(e,s1))),(-6*f,0))
      if not both(s6,*BOUNDS['6'],T,D):ct['s6_reject']+=1;continue
      ct['s6_pass']+=1
      if binary((ONE,A1,A2,A3,A4,A5,A6),(ONE,B1,B2,B3,B4,B5,B6),6,T,N) is None:ct['c6_reject']+=1;continue
      ct['row6_states']+=1
      cross7=ZERO
      Avec=(ONE,A1,A2,A3,A4,A5,A6);Bvec=(ONE,B1,B2,B3,B4,B5,B6)
      for i in range(1,7):cross7=add(cross7,mul(Avec[i],Bvec[7-i],T,N))
      if cross7[1]!=0 or cross7[0] not in (0,1):ct['c7_reject']+=1;continue
      ct['cross7_pass']+=1
      for g in range(-9,10):
       ct['g_box']+=1
       s7=add(add(add(add(add(add(neg(mul(a,s6,T,N)),neg(scale(b,s5))),neg(scale(c,s4))),neg(scale(d,s3))),neg(scale(e,s2))),neg(scale(f,s1))),(-7*g,0))
       if not both(s7,*BOUNDS['7'],T,D):ct['s7_reject']+=1;continue
       s8=add(add(add(add(add(add(add(neg(mul(a,s7,T,N)),neg(scale(b,s6))),neg(scale(c,s5))),neg(scale(d,s4))),neg(scale(e,s3))),neg(scale(f,s2))),neg(scale(g,s1))),(-8,0))
       if not both(s8,*BOUNDS['8'],T,D):ct['s8_reject']+=1;continue
       G7=-g-(b*G5+c*G4+d*G3+e*G2+f*g1)
       Q7=-g-(b*Q5+c*Q4+d*Q3+e*Q2+f*q1)
       # These definitions make A7=B7=0, using the already proved G6=Q6=0.
       assert add((G7+g+b*G5+c*G4+d*G3+e*G2+f*g1,0),ZERO)==ZERO
       assert add((Q7+g+b*Q5+c*Q4+d*Q3+e*Q2+f*q1,0),ZERO)==ZERO
       ct['row7_states']+=1
       final.append([b,c,d,e,f,g,rho,sigma,eta,zeta,alpha,beta,eta5,zeta5,G7,Q7])
  if time.monotonic()-begin>maxseconds:raise TimeoutError(f'core cap T={T},N={N}')
 summary={'T':T,'N':N,'Delta':D,'counts':dict(ct),'row7_states':len(final),'seconds':time.monotonic()-begin}
 atomic(target,{'status':'complete','binding':binding,'summary':summary,'records':final})
 return summary

def symbolic_check()->None:
 import sympy as sp
 a,b,c,d,e,f,g=sp.symbols('a b c d e f g');G1,G2,G3,G4,G5,G6,G7=sp.symbols('G1:8')
 P=[1,a,b,c,d,e,f,g,1];G=[1,G1,G2,G3,G4,G5,G6,G7]
 A=[sp.expand(sum(P[i]*G[k-i] for i in range(k+1))) for k in range(8)]
 assert sp.expand(A[5]-(a*G4+b*G3+c*G2+d*G1+G5+e))==0
 assert sp.expand(A[6]-(G6+a*G5+b*G4+c*G3+d*G2+e*G1+f))==0
 assert sp.expand(A[7]-(a*G6+b*G5+c*G4+d*G3+e*G2+f*G1+g+G7))==0
 # Newton recurrence independently checked through p8=1.
 coeff=[None,a,b,c,d,e,f,g,1];s={0:sp.Integer(8)}
 for k in range(1,9):s[k]=sp.expand(-sum(coeff[j]*s[k-j] for j in range(1,k))-k*coeff[k])
 assert s[8].has(sp.Integer(-8))

def load_stage4(root:Path)->tuple[list[tuple[str,str]],str]:
 summary=root/'summary.json';s=json.loads(summary.read_text());assert s['status']=='complete' and s['cores']==s['completed']==1140 and s['row4_states']==60134
 fs=sorted((root/'cores').glob('*.json'));assert len(fs)==1140
 h=hashlib.sha256();jobs=[];total=0
 for f in fs:
  raw=f.read_bytes();digest=hashlib.sha256(raw).hexdigest();h.update(f.name.encode()+b'\0'+raw)
  x=json.loads(raw);assert x['status']=='complete';total+=len(x['records']);jobs.append((str(f),digest))
 assert total==60134
 return jobs,hashlib.sha256(summary.read_bytes()+h.digest()).hexdigest()

def main()->None:
 ap=argparse.ArgumentParser();ap.add_argument('--stage4',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=1);ap.add_argument('--core-seconds',type=float,default=600);ap.add_argument('--seconds',type=float,default=14400);ap.add_argument('--core-limit',type=int,default=None);ap.add_argument('--only-core',nargs=2,type=int);ap.add_argument('--symbolic-check',action='store_true');args=ap.parse_args()
 if args.symbolic_check: symbolic_check();print(json.dumps({'status':'symbolic_formulas_passed'}));return
 if not __debug__:raise RuntimeError('refuse python -O')
 assert 1<=args.workers<=4
 if args.out.exists() and not (args.out/'cores').is_dir():raise FileExistsError(args.out)
 args.out.mkdir(parents=True,exist_ok=True);(args.out/'cores').mkdir(exist_ok=True)
 jobs,inputhash=load_stage4(args.stage4)
 if args.only_core is not None:
  t,n=args.only_core;jobs=[z for z in jobs if json.loads(Path(z[0]).read_text())['summary']['T']==t and json.loads(Path(z[0]).read_text())['summary']['N']==n]
  assert len(jobs)==1
 if args.core_limit is not None:jobs=jobs[:args.core_limit]
 codehash=sha(Path(__file__));start=time.monotonic();done=[];ct=Counter();deadline=start+args.seconds
 with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
  fut=[pool.submit(one,src,dig,str(args.out),codehash,args.core_seconds) for src,dig in jobs]
  for f in cf.as_completed(fut):
   z=f.result();done.append((z['T'],z['N']));ct.update(z['counts']);print(json.dumps({'done':len(done),'expected':len(jobs),'T':z['T'],'N':z['N'],'row7':z['row7_states']}),flush=True)
   if time.monotonic()>deadline:raise TimeoutError('global cap')
 status='complete' if len(jobs)==1140 else 'partial_scope'
 out={'status':status,'scope':'Necessary S={1}, degree8 rows5--7 from frozen stage4 only; no root screen or arbitrary-cofactor exclusion.','code_sha256':codehash,'stage4_input_sha256':inputhash,'cores':len(jobs),'completed':len(done),'unknown_cores':len(jobs)-len(done),'counts':dict(ct),'row7_states':ct['row7_states'],'workers':args.workers,'seconds':time.monotonic()-start,'python':sys.version}
 atomic(args.out/'summary.json',out);print(json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
