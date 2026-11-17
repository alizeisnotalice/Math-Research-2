#!/usr/bin/env python3
"""Independent degree-seven single-position replay using integer H*G*Q."""
import argparse,hashlib,itertools,json,math,time,sys
from pathlib import Path

def plus(x,y):return x[0]+y[0],x[1]+y[1]
def times(x,y,T,N):return x[0]*y[0]-N*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+T*x[1]*y[1]
def bounded(x,T,D,num,den):
 u,v=x;left=abs(2*u+v*T)*den;right=2*num-left
 return right>0 and v*v*D*den*den<right*right

def floor(x,T,D):
 p,q=x
 if not q:return p
 s=math.isqrt(q*q*D)
 return (2*p+q*T+(s if q>0 else -s-1))//2

def choices(x,T,D):
 f=floor(x,T,D)
 return (-f,1-f) if x[1]==0 else (-f,)

def unit(x,T,D):return x[0] in (0,1) if x[1]==0 else floor(x,T,D)==0

def base(F,G,k,T,conj=False):
 z=sum(F[i]*G[k-i] for i in range(2,k+1));q=G[k-1]
 return (z+q*T,-q) if conj else (z,q)

def ck(F,G,Q,k,T,N):
 H=[sum(F[i]*F[j-i] for i in range(j+1))+ (T*F[j-1] if j else 0)+(N if j==2 else 0) for j in range(k+1)]
 return sum(G[i]*Q[j]*H[k-i-j] for i in range(k+1) for j in range(k-i+1))

def extend(F,states,k,T,N,D):
 out=[]
 for G,Q in states:
  for g,q in itertools.product(choices(base(F,G,k,T),T,D),choices(base(F,Q,k,T,True),T,D)):
   gg=G+(g,);qq=Q+(q,)
   if ck(F,gg,qq,k,T,N) in (0,1):out.append((gg,qq))
 return out

def inputs():
 return sorted((T,N,T*T-4*N) for T in range(-16,17) for N in range(-72,73) if 0<T*T-4*N<(17-abs(T))**2 and math.isqrt(T*T-4*N)**2!=T*T-4*N)

def one(T,N,D):
 start=time.monotonic();c={'c2':0,'c3':0,'c4':0,'prefix4':0,'c5':0,'prefix6':0};passed5=[];passed6=[]
 a=(0,1);a2=times(a,a,T,N);a3=times(a2,a,T,N);a4=times(a3,a,T,N)
 s1=extend((1,0),[((1,),(1,))],1,T,N,D)
 for c2 in range(-29,30):
  if not bounded(plus(a2,(-2*c2,0)),T,D,55,4):continue
  c['c2']+=1;s2=extend((1,0,c2),s1,2,T,N,D)
  for c3 in range(-52,53):
   p3=plus((-a3[0]-3*c3,-a3[1]),(0,3*c2))
   if not bounded(p3,T,D,203,8):continue
   c['c3']+=1;s3=extend((1,0,c2,c3),s2,3,T,N,D)
   for c4 in range(-52,53):
    p4=plus(plus(a4,(-4*c2*a2[0]+2*c2*c2-4*c4,-4*c2*a2[1])),(0,4*c3))
    if not bounded(p4,T,D,787,16):continue
    c['c4']+=1;F=(1,0,c2,c3,c4);s4=extend(F,s3,4,T,N,D);c['prefix4']+=len(s4)
    for G,Q in s4:
     gg=G+(0,);qq=Q+(0,);c5=-(ck(F+(0,),gg,qq,5,T,N)//2)
     if not -29<=c5<=29:continue
     F5=F+(c5,)
     if not unit(base(F5,gg,5,T),T,D) or not unit(base(F5,qq,5,T,True),T,D):continue
     assert ck(F5,gg,qq,5,T,N) in (0,1)
     c['c5']+=1;passed5.append([c2,c3,c4,c5,list(gg),list(qq)])
     # The endpoint G5=Q5=0 makes row 6 rational; choose c6=0 and cancel.
     F6=F5+(0,);g6=-base(F6,gg+(0,),6,T)[0];q6=-base(F6,qq+(0,),6,T,True)[0]
     gg6=gg+(g6,);qq6=qq+(q6,)
     if ck(F6,gg6,qq6,6,T,N) in (0,1):
      c['prefix6']+=1;passed6.append([c2,c3,c4,c5,list(gg6),list(qq6)])
 return {'TN':[T,N,D],'counts':c,'passed5':passed5,'passed6':passed6,'elapsed_seconds':time.monotonic()-start}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--T',type=int,required=True);ap.add_argument('--N',type=int,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();triple=(a.T,a.N,a.T*a.T-4*a.N);assert triple in inputs();r=one(*triple);r.update(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=sys.version);a.out.write_text(json.dumps(r,sort_keys=True)+'\n');print(json.dumps({k:r[k] for k in ['TN','counts','elapsed_seconds']}))
if __name__=='__main__':main()
