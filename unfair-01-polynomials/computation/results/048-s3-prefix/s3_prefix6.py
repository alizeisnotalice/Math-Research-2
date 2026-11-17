#!/Users/zzh/.venvs/main/bin/python
import json,hashlib
from math import isqrt
from pathlib import Path
H=Path(__file__).resolve().parent
def add(*z):return(sum(x[0] for x in z),sum(x[1] for x in z))
def mul(x,y,T,N):return(x[0]*y[0]-N*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+T*x[1]*y[1])
def fl(z,T,D):
 p,q=z
 if not q:return p
 u=2*p+q*T;r=isqrt(q*q*D);return (u+r)//2 if q>0 else (u-r-1)//2
def choices(z,T,D):
 return (-z[0],1-z[0]) if z[1]==0 else (-fl(z,T,D),)
def unit(z,T,D):return z[0] in (0,1) if z[1]==0 else fl(z,T,D)==0
def coeff(P,G,k,T,N):return add(*(mul(P[j],(G[k-j],0),T,N) for j in range(min(k,6)+1)))
def C(A,B,k,T,N):z=add(*(mul(A[i],B[k-i],T,N) for i in range(k+1)));return z[1]==0 and z[0] in (0,1)
def tn():
 for T in range(-16,17):
  R=203-12*abs(T)
  for D in range(1,(R*R-1)//144+1):
   if (T*T-D)%4==0 and isqrt(D)**2!=D:yield T,(T*T-D)//4,D
out=[]; counts={'TN':0,'cuv':0,'C4':0,'C5':0,'C6':0}
for T,N,D in tn():
 counts['TN']+=1
 for c in range(-6,7):
  for u,v in ((0,0),(1,0),(0,1)):
   counts['cuv']+=1; ga=-fl((0,1),T,D); gb=-fl((T,-1),T,D)
   G=(1,0,u-c,ga,-c*(u-c));Q=(1,0,v-c,gb,-c*(v-c))
   P=((1,0),(0,0),(c,0),(0,1),(0,0),(0,0),(0,0));Pb=((1,0),(0,0),(c,0),(T,-1),(0,0),(0,0),(0,0))
   A=[coeff(P,G,k,T,N) for k in range(5)];B=[coeff(Pb,Q,k,T,N) for k in range(5)]
   if not(all(unit(z,T,D) for z in A+B) and all(C(A,B,k,T,N) for k in range(5))):continue
   counts['C4']+=1
   for g5 in choices(coeff(P,G+(0,),5,T,N),T,D):
    for q5 in choices(coeff(Pb,Q+(0,),5,T,N),T,D):
     G5=G+(g5,);Q5=Q+(q5,);A=[coeff(P,G5,k,T,N) for k in range(6)];B=[coeff(Pb,Q5,k,T,N) for k in range(6)]
     if not(all(unit(z,T,D) for z in A+B) and all(C(A,B,k,T,N) for k in range(6))):continue
     counts['C5']+=1
     for g6 in choices(coeff(P,G5+(0,),6,T,N),T,D):
      for q6 in choices(coeff(Pb,Q5+(0,),6,T,N),T,D):
       GG=G5+(g6,);QQ=Q5+(q6,);A=[coeff(P,GG,k,T,N) for k in range(7)];B=[coeff(Pb,QQ,k,T,N) for k in range(7)]
       if not(all(unit(z,T,D) for z in A+B) and all(C(A,B,k,T,N) for k in range(7))):continue
       counts['C6']+=1;out.append([T,N,D,c,u,v,*GG,*QQ])
r={'status':'PREFIX6_COMPLETE_REPRESENTATIVE','counts':counts,'states':out,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(H/'s3_prefix6.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'counts':counts,'states':len(out)}))
