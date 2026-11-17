#!/Users/zzh/.venvs/main/bin/python
"""Exact Newton necessary screen for degree-7 S={3}; shardable, no prefix search."""
import argparse,json,hashlib
from math import isqrt
from pathlib import Path
H=Path(__file__).resolve().parent
def add(*z):return(sum(x[0] for x in z),sum(x[1] for x in z))
def mul(x,y,T,N):return(x[0]*y[0]-N*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+T*x[1]*y[1])
def bound(z,T,D,num,den):
 p,q=z; R=2*num-den*abs(2*p+q*T); return R>0 and (den*q)**2*D<R*R
def tn():
 for T in range(-16,17):
  R=203-12*abs(T)
  for D in range(1,(R*R-1)//144+1):
   if (T*T-D)%4==0 and isqrt(D)**2!=D:yield T,(T*T-D)//4,D
def run(start,count,keep=True):
 C={'TN':0,'c':0,'d':0,'e':0,'s7':0,'f':0};out=[]
 for T,N,D in list(tn())[start:start+count]:
  C['TN']+=1
  for c in range(-6,7):
   C['c']+=1
   for d in range(-52,53):
    if abs(2*c*c-4*d)*16>=787:continue
    C['d']+=1
    for e in range(-29,30):
     if not bound((-5*e,5*c),T,D,3107,32):continue
     C['e']+=1
     if not bound((7*(c*e-1),7*(d-c*c)),T,D,49283,128):continue
     C['s7']+=1
     # s6=3a^2-2c^3+6cd-6f.
     s6=add(mul((3,0),( -N,T),T,N),(-2*c**3+6*c*d-6*f,0)) if False else None
     for f in range(-8,9):
      z=add(mul((3,0),(-N,T),T,N),(-2*c**3+6*c*d-6*f,0))
      if bound(z,T,D,12355,64):
       C['f']+=1
       if keep:out.append([T,N,D,c,d,e,f])
 return C,out
def main():
 p=argparse.ArgumentParser();p.add_argument('--tn-start',type=int,default=0);p.add_argument('--tn-count',type=int,default=1);p.add_argument('--counts-only',action='store_true');z=p.parse_args();C,out=run(z.tn_start,z.tn_count,not z.counts_only);r={'status':'PARTIAL_NEWTON_SCREEN','start':z.tn_start,'count':z.tn_count,'TN_total':len(list(tn())),'counts':C,'mothers':out,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};q=H/f's3_newton_{z.tn_start}_{z.tn_count}.json';q.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'counts':C,'mothers':len(out)}))
if __name__=='__main__':main()
