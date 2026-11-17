#!/Users/zzh/.venvs/main/bin/python
"""Exact low-prefix necessary screen for degree-7 S={1}; batchable by T index."""
import argparse,json,hashlib
from math import isqrt
from pathlib import Path
H=Path(__file__).resolve().parent
def fl(p,q,T,D):
 if not q:return p
 u=2*p+q*T;r=isqrt(q*q*D)
 return (u+r)//2 if q>0 else (u-r-1)//2
def opt(base,T,D):
 p,q=base;z=fl(p,q,T,D);return (-z,1-z) if q==0 else (-z,)
def add(*z):return (sum(x[0] for x in z),sum(x[1] for x in z))
def mul(x,y,T,N):return (x[0]*y[0]-N*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+T*x[1]*y[1])
def bounded(x,T,D,num,den):
 # both |p+q a|<num/den
 p,q=x; A=den*abs(2*p+q*T); B=den*abs(q); R=2*num-A
 return R>0 and B*B*D<R*R
def coeff(P,G,k,T,N):return add(*(mul(P[j],(G[k-j],0),T,N) for j in range(min(k,len(P)-1)+1)))
def base(P,G,k):
 return (sum(P[j][0]*G[k-j] for j in range(1,min(k,len(P)-1)+1)),sum(P[j][1]*G[k-j] for j in range(1,min(k,len(P)-1)+1)))
def validC(A,B,k,T,N):
 z=add(*(mul(A[i],B[k-i],T,N) for i in range(k+1)));return z[1]==0 and z[0] in (0,1)
def unit(z,T,D):
 return z[0] in (0,1) if z[1]==0 else fl(*z,T,D)==0
def tn():
 for T in range(-16,17):
  for D in range(1,(17-abs(T))**2):
   if (T*T-D)%4==0 and isqrt(D)**2!=D:yield T,(T*T-D)//4,D
def run(start,count):
 counts={'TN':0,'c2':0,'c3':0,'c4':0,'prefix4':0,'c5':0,'prefix6':0}; survivors=[]
 rows=list(tn())[start:start+count]
 for ix,(T,N,D) in enumerate(rows):
  counts['TN']+=1
  a2=(-N,T);a3=(-N*T,T*T-N);a4=(N*N-N*T*T,T**3-2*N*T)
  for c2 in range(-29,30):
   if not bounded(add(a2,(-2*c2,0)),T,D,55,4):continue
   counts['c2']+=1
   for c3 in range(-52,53):
    if not bounded(add((-a3[0]-3*c3,a3[1]*-1), (0,3*c2)),T,D,203,8):continue
    counts['c3']+=1
    for c4 in range(-52,53):
     p4=add(a4,(-4*c2*a2[0]+2*c2*c2-4*c4,-4*c2*a2[1]),(0,4*c3))
     if not bounded(p4,T,D,787,16):continue
     counts['c4']+=1
     # P and conjugate in the common basis a; c5 filled below.
     P0=((1,0),(0,1),(c2,0),(c3,0),(c4,0),(0,0),(0,0)); Q0=((1,0),(T,-1),(c2,0),(c3,0),(c4,0),(0,0),(0,0)); states=[((1,),(1,))]
     for k in range(1,5):
      nxt=[]
      for G,Q in states:
       for gi in opt(base(P0,G,k),T,D):
        for qi in opt(base(Q0,Q,k),T,D):
         gg=G+(gi,);qq=Q+(qi,);A=[coeff(P0,gg,i,T,N) for i in range(k+1)];B=[coeff(Q0,qq,i,T,N) for i in range(k+1)]
         if all(validC(A,B,i,T,N) for i in range(k+1)):nxt.append((gg,qq))
      states=nxt
     counts['prefix4']+=len(states)
     for G,Q in states:
      c5set=set(opt(base(P0,G,5),T,D))&set(opt(base(Q0,Q,5),T,D))&set(range(-29,30))
      for c5 in c5set:
       P=P0[:5]+((c5,0),(0,0));QP=Q0[:5]+((c5,0),(0,0));gg=G+(0,);qq=Q+(0,)
       A=[coeff(P,gg,i,T,N) for i in range(6)];B=[coeff(QP,qq,i,T,N) for i in range(6)]
       if not all(unit(z,T,D) for z in A+B) or not all(validC(A,B,i,T,N) for i in range(6)):continue
       counts['c5']+=1
       # Set c6=0 temporarily: forced G6,Q6 cancel it, while C6 is unchanged.
       g6=-coeff(P,gg+(0,),6,T,N)[0];q6=-coeff(QP,qq+(0,),6,T,N)[0];gg=gg+(g6,);qq=qq+(q6,)
       A=[coeff(P,gg,i,T,N) for i in range(7)];B=[coeff(QP,qq,i,T,N) for i in range(7)]
       if all(unit(z,T,D) for z in A+B) and all(validC(A,B,i,T,N) for i in range(7)):
        counts['prefix6']+=1;survivors.append({'T':T,'N':N,'D':D,'c2':c2,'c3':c3,'c4':c4,'c5':c5,'c6':'unfixed[-8,8]','G':list(gg),'Q':list(qq)})
 return counts,survivors
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--tn-start',type=int,default=0);ap.add_argument('--tn-count',type=int,default=1);z=ap.parse_args();c,s=run(z.tn_start,z.tn_count);o={'status':'PARTIAL_TN_PREFIX_SCREEN','tn_start':z.tn_start,'tn_count':z.tn_count,'input_TN_total':len(list(tn())),'counts':c,'survivors':s,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};p=H/f's1_shard_{z.tn_start}_{z.tn_count}.json';p.write_text(json.dumps(o,indent=2)+'\n');print(json.dumps({'counts':c,'survivors':len(s)}))
if __name__=='__main__':main()
