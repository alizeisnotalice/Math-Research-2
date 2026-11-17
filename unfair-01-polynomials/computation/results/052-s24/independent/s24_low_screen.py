#!/usr/bin/env python3
"""Exact low-row sieve for degree-seven S={2,4}; no e/f scan.

P=1+a*x^2+c*x^3+h*x^4+e*x^5+f*x^6+x^7.  The nonrational a,h
are independent full-O_K elements.  This enumerates all lattice pairs after
the forced s2 bound, filters s4, and explicitly records every c=-8..8 row
surviving through C4=1.  It assumes no cofactor degree bound.
"""
import argparse,hashlib,json,math,os,sys,time
from collections import Counter,defaultdict
from pathlib import Path

def sf(d):
 p=2
 while p*p<=d:
  if d%(p*p)==0:return False
  p+=1
 return True
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def scl(k,x):return(k*x[0],k*x[1])
def mul(x,y,D):
 m,n=x;u,v=y
 return (m*u+n*v*((D-1)//4),m*v+n*u+n*v) if D%4==1 else (m*u+D*n*v,m*v+n*u)
def bar(x,D):return(x[0]+x[1],-x[1]) if D%4==1 else(x[0],-x[1])
def rv(x,D):return(2*x[0]+x[1],x[1],2) if D%4==1 else(x[0],x[1],1)
def sign(x,D):
 R,V,s=rv(x,D)
 if V==0:return(R>0)-(R<0)
 if V>0:
  if R>=0:return 1
  return 1 if V*V*D>R*R else -1
 if R<=0:return-1
 return 1 if R*R>V*V*D else -1
def abslim(x,D,num,den):
 R,V,s=rv(x,D);room=s*num-den*abs(R)
 return room>0 and (den*V)*(den*V)*D<room*room
def floorq(x,D):
 lo,hi=-1,1
 while sign(add(x,(-lo,0)),D)<0:hi=lo;lo*=2
 while sign(add(x,(-hi,0)),D)>=0:lo=hi;hi*=2
 while hi-lo>1:
  mid=(lo+hi)//2
  if sign(add(x,(-mid,0)),D)>=0:lo=mid
  else:hi=mid
 return lo
def one_shift(base,upper,D):
 lo=-floorq(base,D);hi=floorq(add(upper,neg(base)),D)
 if lo>hi:return[]
 if lo!=hi:raise AssertionError(('nonunique interval',base,upper,lo,hi))
 z=lo;v=add(base,(z,0));assert sign(v,D)>=0 and sign(add(upper,neg(v)),D)>=0
 return[z]
def tr(x,D):return 2*x[0]+x[1] if D%4==1 else 2*x[0]
def cand(D,num,den,positive):
 lim=(2*num+den-1)//den if D%4==1 else (num+den-1)//den
 z=[]
 for n in range(1 if positive else -lim,lim+1):
  if n==0:continue
  for m in range(-lim,lim+1):
   if abslim((m,n),D,num,den):z.append((m,n))
 return z
def fields():
 out=[]
 # s2=-2a makes this range complete; keeping 189 is a harmless exact cap.
 for D in range(2,190):
  if not sf(D):continue
  aa=cand(D,55,8,True)
  if aa:out.append((D,aa,cand(D,419,8,False)))
 return out
def atom(path,data):
 t=path.with_suffix(path.suffix+'.partial')
 with t.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(t,path)

def screen_pair(D,a,h):
 """Return base low row or a rejection tag. h is p4, not field discriminant."""
 # p1=0: s2=-2a is already enforced in field construction.
 s4=add(scl(2,mul(a,a,D)),scl(-4,h))
 if not abslim(s4,D,787,16):return None,'s4'
 T=tr(a,D);g2=-floorq(a,D);q2=1-T-g2
 A2=add(a,(g2,0));B2=add(bar(a,D),(q2,0))
 if not(add(A2,B2)==(1,0) and sign(A2,D)>0 and sign(B2,D)>0):raise AssertionError(('row2',D,a,h,g2,q2,A2,B2))
 S4=add((1,0),neg(mul(A2,B2,D)))
 if not(sign(S4,D)>0 and sign(add(S4,(-1,0)),D)<0):raise AssertionError(('S4',D,a,h,S4))
 # G1=Q1=0, G3=Q3=-c.  Thus c does not occur in row four.
 baseA=add(h,scl(g2,a));baseB=add(bar(h,D),scl(q2,bar(a,D)))
 shifts=one_shift(baseA,S4,D)
 if not shifts:return None,'A4'
 g4=shifts[0]
 q4v=add(add(S4,neg(baseA)),neg(baseB));q4v=add(q4v,(-g4,0))
 if q4v[1]!=0:raise AssertionError(('q4 irrational',D,a,h,q4v))
 q4=q4v[0]
 A4=add(baseA,(g4,0));B4=add(baseB,(q4,0))
 if not(add(A4,B4)==S4 and sign(A4,D)>=0 and sign(B4,D)>=0):raise AssertionError(('row4',D,a,h,g4,q4,A4,B4,S4))
 return {'D':D,'a_m':a[0],'a_n':a[1],'h_m':h[0],'h_n':h[1],'T':T,
         'g0':1,'g1':0,'g2':g2,'q0':1,'q1':0,'q2':q2,'g4':g4,'q4':q4},'base'

def one_field(D,aa,hh):
 ct=Counter(raw_ad=0);rows=[]
 for a in aa:
  for h in hh:
   ct['raw_ad']+=1;r,code=screen_pair(D,a,h);ct[code]+=1
   if r is not None:
    # All 17 c values are retained exactly as an affine family rather than
    # materialised as 17 identical base records.  A3=B3=0 iff g3=q3=-c,
    # and row four is c-independent in this support type.
    x=dict(r);x.update({'c_min':-8,'c_max':8,'g3_equals':'-c','q3_equals':'-c','C2':1,'C3':0,'C4':1});rows.append(x);ct['c_lifts']+=17
 assert ct['raw_ad']==len(aa)*len(hh)
 return dict(ct),rows

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--start',type=int,default=0);ap.add_argument('--count',type=int,default=None);ap.add_argument('--seconds',type=float,default=14400);z=ap.parse_args();z.out.mkdir(parents=True,exist_ok=True)
 fs=fields();assert sum(len(a)*len(h) for _,a,h in fs)==774100
 if z.count is None:z.count=len(fs)
 assert 0<=z.start and 0<=z.count and z.start+z.count<=len(fs)
 source=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();started=time.monotonic();done=[]
 for i in range(z.start,z.start+z.count):
  D,a,h=fs[i];path=z.out/f'field-{i:03d}-D{D}.json';bind={'source_sha256':source,'index':i,'D':D,'a_count':len(a),'h_count':len(h)}
  if path.exists():
   old=json.loads(path.read_text());assert old['status']=='complete' and old['binding']==bind;done.append(i);continue
  if time.monotonic()-started>z.seconds:break
  then=time.monotonic();ct,rows=one_field(D,a,h);obj={'status':'complete','binding':bind,'counts':ct,'rows':rows,'elapsed_seconds':time.monotonic()-then};atom(path,(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode());done.append(i);print(json.dumps({'index':i,'D':D,'counts':ct,'rows':len(rows),'seconds':obj['elapsed_seconds']},sort_keys=True),flush=True)
 total=Counter();rows=[]
 for i in done:
  D=fs[i][0];x=json.loads((z.out/f'field-{i:03d}-D{D}.json').read_text());total.update(x['counts']);rows.extend(x['rows'])
 out={'status':'complete' if len(done)==z.count else 'partial_unknown','scope':'S={2,4} exact rows through C4 only; c expanded, e/f unenumerated and no cofactor degree bound.','completed_indices':done,'counts':dict(total),'rows':rows,'source_sha256':source,'python':sys.version,'elapsed_seconds':time.monotonic()-started}
 atom(z.out/'summary.json',(json.dumps(out,sort_keys=True,indent=2)+'\n').encode());print(json.dumps({'status':out['status'],'fields':len(done),'counts':dict(total),'rows':len(rows),'seconds':out['elapsed_seconds']},sort_keys=True))
if __name__=='__main__':main()
