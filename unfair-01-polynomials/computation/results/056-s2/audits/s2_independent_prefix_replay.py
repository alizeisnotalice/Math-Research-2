#!/usr/bin/env python3
"""Independent exact S={2}, degree-8 necessary-prefix replay.

Representation is (u+v sqrt(D))/2, deliberately different from production's
p+q*a.  It rebuilds every A/B/C row directly and retains rational endpoints.
Only necessary prefix conditions are tested.
"""
from __future__ import annotations
import argparse, concurrent.futures as cf, hashlib, json, math, os, time
from collections import Counter
from pathlib import Path

Z=(0,0); O=(2,0)
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def neg(x): return (-x[0],-x[1])
def sc(n,x): return (n*x[0],n*x[1])
def ci(n): return (2*n,0)
def mul(x,y,D):
 u,v=x;r,s=y; a=u*r+v*s*D; b=u*s+v*r
 assert a%2==b%2==0
 return (a//2,b//2)
def conj(x): return (x[0],-x[1])
def sgnraw(u,v,D):
 if v==0:return (u>0)-(u<0)
 if u==0:return (v>0)-(v<0)
 if u*v>0:return (u>0)-(u<0)
 q=(u*u>v*v*D)-(u*u<v*v*D)
 return q if u>0 else -q
def sign(x,D): return sgnraw(x[0],x[1],D)
def lt(x,num,den,D): return sgnraw(den*x[0]-2*num,den*x[1],D)<0
def bothabs(x,num,den,D): return lt(x,num,den,D) and lt(neg(x),num,den,D) and lt(conj(x),num,den,D) and lt(neg(conj(x)),num,den,D)
def floorq(x,D):
 R=abs(x[0])+abs(x[1])*(math.isqrt(D)+3)+8; lo,hi=-R,R
 while hi-lo>1:
  m=(lo+hi)//2
  if sign((x[0]-2*m,x[1]),D)>=0:lo=m
  else:hi=m
 return lo
def shifts(x,D):
 if x[1]==0:
  assert x[0]%2==0
  n=x[0]//2;return [(-n,Z),(1-n,O)]
 h=-floorq(x,D); y=(x[0]+2*h,x[1])
 assert sign(y,D)>0 and sign(add(O,neg(y)),D)>0
 return [(h,y)]
def unit(x): return x[1]==0 and x[0] in (0,2)
def convcoef(A,B,k,D):
 z=Z
 for i in range(k+1):z=add(z,mul(A[i],B[k-i],D))
 return z
def newton(P,D,through):
 s=[ci(8)]
 for k in range(1,through+1):
  z=Z
  for j in range(1,k):z=add(z,mul(P[j],s[k-j],D))
  s.append(neg(add(z,sc(k,P[k]))))
 return s
def cores():
 out=[]
 for L in range(-9,10):
  for T in range(-82,83):
   x=17-abs(L*L-T);y=83-abs(T)
   if x<=0 or y<=0:continue
   for D in range(1,min(x*x,y*y)):
    if math.isqrt(D)**2==D or D%4!=(T*T)%4:continue
    N=(T*T-D)//4
    if bothabs((T,1),83,2,D):out.append((L,T,N,D))
 return sorted(out)
def atomic(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix('.partial')
 with t.open('w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def one4(job):
 L,T,N,D,out,cap=job;ct=Counter(); rec=[];beg=time.monotonic(); a=(T,1)
 p=Path(out)/'stage4'/f'L{L}_T{T}_N{N}.json'
 for A1n,B1n in ((0,0),(1,0),(0,1)):
  A1=ci(A1n);B1=ci(B1n);g1=A1n-L;q1=B1n-L
  for g2,A2 in shifts(add(a,ci(L*g1)),D):
   for q2,B2 in shifts(add(conj(a),ci(L*q1)),D):
    if not unit(convcoef([O,A1,A2],[O,B1,B2],2,D)):continue
    ct['row2']+=1
    for b in range(-92,93):
     P=[O,ci(L),a,ci(b),Z,Z,Z,Z,O];ss=newton(P,D,3)
     if not bothabs(ss[3],65,2,D):continue
     baseA=add(ci(b+L*g2),sc(g1,a));baseB=add(ci(b+L*q2),sc(q1,conj(a)))
     for g3,A3 in shifts(baseA,D):
      for q3,B3 in shifts(baseB,D):
       if not unit(convcoef([O,A1,A2,A3],[O,B1,B2,B3],3,D)):continue
       ct['row3']+=1
       for c in range(-120,121):
        P[4]=ci(c);ss=newton(P,D,4)
        if not bothabs(ss[4],257,4,D):continue
        A4=add(ci(c+L*g3+b*g1),add(sc(g2,a),Z));B4=add(ci(c+L*q3+b*q1),sc(q2,conj(a)))
        # High endpoint has already forced G4=Q4=0.
        if not (sign(A4,D)>=0 and sign(add(O,neg(A4)),D)>=0 and sign(B4,D)>=0 and sign(add(O,neg(B4)),D)>=0):continue
        if not unit(convcoef([O,A1,A2,A3,A4],[O,B1,B2,B3,B4],4,D)):continue
        ct['row4']+=1;rec.append([L,T,N,D,b,c,g1,q1,g2,q2,g3,q3])
  if time.monotonic()-beg>cap:raise TimeoutError((L,T,N,D))
 atomic(p,{'status':'complete','core':[L,T,N,D],'counts':dict(ct),'records':rec});return dict(ct),len(rec)
def one58(job):
 src,out,cap=job;x=json.loads(Path(src).read_text());L,T,N,D=x['core'];a=(T,1);ct=Counter();rec=[];beg=time.monotonic();p=Path(out)/'stage58'/Path(src).name
 for L,T,N,D,b,c,g1,q1,g2,q2,g3,q3 in x['records']:
  A1=ci(L+g1);B1=ci(L+q1);A2=add(a,ci(L*g1+g2));B2=add(conj(a),ci(L*q1+q2));A3=add(ci(b+L*g2+g3),sc(g1,a));B3=add(ci(b+L*q2+q3),sc(q1,conj(a)));A4=add(ci(c+L*g3+b*g1),sc(g2,a));B4=add(ci(c+L*q3+b*q1),sc(q2,conj(a)))
  for d in range(-92,93):
   P=[O,ci(L),a,ci(b),ci(c),ci(d),Z,Z,O];ss=newton(P,D,5)
   if not bothabs(ss[5],1025,8,D):continue
   baseA=add(ci(d+b*g2+c*g1),sc(g3,a));baseB=add(ci(d+b*q2+c*q1),sc(q3,conj(a)))
   for g5,A5 in shifts(baseA,D):
    for q5,B5 in shifts(baseB,D):
     if not unit(convcoef([O,A1,A2,A3,A4,A5],[O,B1,B2,B3,B4,B5],5,D)):continue
     ct['row5']+=1
     cross6=convcoef([O,A1,A2,A3,A4,A5,Z],[O,B1,B2,B3,B4,B5,Z],6,D)
     if not unit(cross6):continue
     for e in range(-41,42):
      P[6]=ci(e);ss=newton(P,D,6)
      if not bothabs(ss[6],4097,16,D):continue
      vA=ci(e+L*g5+b*g3+c*g2+d*g1);vB=ci(e+L*q5+b*q3+c*q2+d*q1)
      # aG4 and aQ4 vanish; A6=B6=0 fixes G6,Q6 separately.
      if vA[1] or vB[1] or vA[0]%2 or vB[0]%2:continue
      g6=-vA[0]//2;q6=-vB[0]//2;ct['row6']+=1
      for f in range(-9,10):
       P[7]=ci(f);ss=newton(P,D,8)
       if not bothabs(ss[7],16385,32,D) or not bothabs(ss[8],65537,64,D):continue
       baseA=add(ci(f+L*g6+b*0+c*g3+d*g2+e*g1),sc(g5,a));baseB=add(ci(f+L*q6+c*q3+d*q2+e*q1),sc(q5,conj(a)))
       for g7,A7 in shifts(baseA,D):
        for q7,B7 in shifts(baseB,D):
         if not unit(convcoef([O,A1,A2,A3,A4,A5,Z,A7],[O,B1,B2,B3,B4,B5,Z,B7],7,D)):continue
         ct['row7']+=1;rec.append([L,T,N,D,b,c,d,e,f,g1,q1,g2,q2,g3,q3,g5,q5,g6,q6,g7,q7])
  if time.monotonic()-beg>cap:raise TimeoutError((L,T,N,D))
 atomic(p,{'status':'complete','core':[L,T,N,D],'counts':dict(ct),'records':rec});return dict(ct),len(rec)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=('4','58'));ap.add_argument('--out',required=True);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--seconds',type=float,default=14400);a=ap.parse_args();out=Path(a.out)
 if a.stage=='4': jobs=[(*z,str(out),a.seconds) for z in cores()]
 else: jobs=[(str(p),str(out),a.seconds) for p in sorted((out/'stage4').glob('*.json'))]
 assert len(jobs)==11760
 C=Counter();n=0
 with cf.ProcessPoolExecutor(max_workers=a.workers) as ex:
  for d,m in ex.map(one4 if a.stage=='4' else one58,jobs):C.update(d);n+=m
 atomic(out/f'summary-{a.stage}.json',{'status':'complete','stage':a.stage,'cores':len(jobs),'states':n,'counts':dict(C)})
 print(json.dumps({'status':'complete','stage':a.stage,'states':n,'counts':dict(C)},sort_keys=True))
if __name__=='__main__':main()
