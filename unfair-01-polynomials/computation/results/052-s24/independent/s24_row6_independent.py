#!/usr/bin/env python3
"""Independent sixth-row shift sieve for S={2,4} base records."""
import argparse,hashlib,json,math,os,sys,time
from collections import Counter,defaultdict
from pathlib import Path
INPUT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/round052-s24/low-full/summary.json')
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def scl(k,x):return(k*x[0],k*x[1])
def mul(x,y,D):
 m,n=x;u,v=y
 return(m*u+n*v*((D-1)//4),m*v+n*u+n*v) if D%4==1 else(m*u+D*n*v,m*v+n*u)
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
 if lo!=hi:raise AssertionError(('nonunique',base,upper,lo,hi))
 return[lo]
def atom(path,data):
 t=path.with_suffix(path.suffix+'.partial')
 with t.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(t,path)
def one(r):
 D=r['D'];a=(r['a_m'],r['a_n']);h=(r['h_m'],r['h_n']);aa,hh=bar(a,D),bar(h,D)
 g2,q2,g4,q4=r['g2'],r['q2'],r['g4'],r['q4']
 A2=add(a,(g2,0));B2=add(aa,(q2,0));A4=add(add(h,scl(g2,a)),(g4,0));B4=add(add(hh,scl(q2,aa)),(q4,0))
 cross=add(mul(A2,B4,D),mul(A4,B2,D));target=add((1,0),neg(cross))
 if not(sign(cross,D)>0 and sign(target,D)>0 and sign(add(target,(-1,0)),D)<0):raise AssertionError(('cross',r,cross,target))
 baseA=add(scl(g4,a),scl(g2,h));baseB=add(scl(q4,aa),scl(q2,hh))
 z=one_shift(baseA,target,D)
 if not z:return None,'R6'
 R6=z[0];S6v=add(add(target,neg(baseA)),neg(baseB));S6v=add(S6v,(-R6,0))
 if S6v[1]!=0:raise AssertionError(('S6 irrational',r,S6v))
 S6=S6v[0];A6=add(baseA,(R6,0));B6=add(baseB,(S6,0))
 if not(add(A6,B6)==target and sign(A6,D)>=0 and sign(B6,D)>=0):raise AssertionError(('row6',r,R6,S6,A6,B6,target))
 out=dict(r);out.update({'R6_equals':'G6+f-c^2','S6_equals':'Q6+f-c^2','R6':R6,'S6':S6,'G6_equals':'R6-f+c^2','Q6_equals':'S6-f+c^2','C6':1})
 return out,'base'
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--start',type=int,default=0);ap.add_argument('--count',type=int,default=None);ap.add_argument('--seconds',type=float,default=14400);z=ap.parse_args();z.out.mkdir(parents=True,exist_ok=True)
 ibytes=INPUT.read_bytes();inp=json.loads(ibytes);assert inp['status']=='complete' and len(inp['rows'])==33575
 groups=defaultdict(list)
 for r in inp['rows']:groups[r['D']].append(r)
 fs=sorted(groups);assert len(fs)==53
 if z.count is None:z.count=len(fs)
 assert 0<=z.start and 0<=z.count and z.start+z.count<=len(fs)
 source=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();bind0={'source_sha256':source,'input_sha256':hashlib.sha256(ibytes).hexdigest()};started=time.monotonic();done=[]
 for i in range(z.start,z.start+z.count):
  D=fs[i];path=z.out/f'field-{i:03d}-D{D}.json';bind={**bind0,'index':i,'D':D,'bases':len(groups[D])}
  if path.exists():
   x=json.loads(path.read_text());assert x['status']=='complete' and x['binding']==bind;done.append(i);continue
  if time.monotonic()-started>z.seconds:break
  ct=Counter();rows=[];then=time.monotonic()
  for r in groups[D]:
   o,c=one(r);ct[c]+=1
   if o is not None:rows.append(o)
  obj={'status':'complete','binding':bind,'counts':dict(ct),'rows':rows,'elapsed_seconds':time.monotonic()-then};atom(path,(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode());done.append(i);print(json.dumps({'index':i,'D':D,'counts':dict(ct),'rows':len(rows),'seconds':obj['elapsed_seconds']},sort_keys=True),flush=True)
 total=Counter();rows=[]
 for i in done:
  D=fs[i];x=json.loads((z.out/f'field-{i:03d}-D{D}.json').read_text());total.update(x['counts']);rows.extend(x['rows'])
 out={'status':'complete' if len(done)==z.count else 'partial_unknown','scope':'S={2,4} exact row6 shift only; c/f affine parameters retained, e unenumerated.','completed_indices':done,'counts':dict(total),'rows':rows,'binding':bind0,'elapsed_seconds':time.monotonic()-started,'python':sys.version}
 atom(z.out/'summary.json',(json.dumps(out,sort_keys=True,indent=2)+'\n').encode());print(json.dumps({'status':out['status'],'fields':len(done),'counts':dict(total),'rows':len(rows),'seconds':out['elapsed_seconds']},sort_keys=True))
if __name__=='__main__':main()
