#!/usr/bin/env python3
"""Exact necessary root-annulus screen for stage-5--8 S={2} mothers.

For P=F+a x^2, this independently builds H=P*Pbar=F^2+T x^2 F+N x^4
in Z[x], then checks positive-axis samples and strict outer/inner Schur
descent.  A rejection is necessary only; a pass is never a certificate.
"""
from __future__ import annotations
import argparse,concurrent.futures as cf,hashlib,json,math,os,sys,time
from collections import Counter
from pathlib import Path
POINTS=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.partial')
 with t.open('w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def conv(a,b):
 z=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i+j]+=x*y
 return z
def H(row):
 L,T,N,D,b,c,d,e,f=row;F=[1,L,0,b,c,d,e,f,1];h=conv(F,F)
 for i,x in enumerate(F):h[i+2]+=T*x
 h[4]+=N
 assert len(h)==17 and h[0]==h[-1]==1
 return h
def schur(a):
 """first failed degree, or None precisely on strict unit disk."""
 a=list(a)
 while len(a)>1:
  g=math.gcd(*a)
  if a[0]<0:g=-g
  a=[x//g for x in a]
  if abs(a[-1])>=a[0]:return len(a)-1
  lead,const=a[0],a[-1];n=len(a)-1
  a=[lead*a[i]-const*a[n-i] for i in range(n)]
 return None
def screen(row):
 h=H(row)
 for j,(p,q) in enumerate(POINTS):
  if sum(v*p**i*q**(16-i) for i,v in enumerate(h))<=0:return 10+j,h
 k=schur([v*(2**i) for i,v in enumerate(h)][::-1])
 if k is not None:return 30+k,h
 k=schur([v*(2**(16-i)) for i,v in enumerate(h)])
 if k is not None:return 50+k,h
 return 255,h
def one(job):
 source,source_hash,outdir,input_hash,codehash,maxseconds=job;x=json.loads(Path(source).read_text());core=x['binding']['core'];started=time.monotonic();target=Path(outdir)/'cores'/Path(source).name
 bind={'version':1,'core':core,'stage58_source_sha256':source_hash,'stage58_input_sha256':input_hash,'code_sha256':codehash,'points':[list(q) for q in POINTS]}
 if target.exists():
  o=json.loads(target.read_text())
  if o.get('binding')==bind and o.get('status')=='complete':return o['summary']
  raise RuntimeError('stale root cache '+str(target))
 rows={tuple(r[:9]) for r in x['records']};ct=Counter();out=[]
 for row in sorted(rows):
  code,h=screen(row);ct['passed' if code==255 else 'rejected']+=1;ct[f'code_{code}']+=1
  out.append({'row':list(row),'H':h,'code':code,'status':'passed' if code==255 else 'rejected'})
  if time.monotonic()-started>maxseconds:raise TimeoutError('root cap '+str(core))
 summary={'core':core,'mothers':len(rows),'counts':dict(ct),'seconds':time.monotonic()-started}
 atomic(target,{'status':'complete','binding':bind,'summary':summary,'records':out});return summary
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stage58',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=1);ap.add_argument('--max-seconds',type=float,default=14400);ap.add_argument('--only-core');a=ap.parse_args()
 s=json.loads((a.stage58/'summary.json').read_text());assert s['status']=='complete' and s['cores']==s['completed']==11760
 fs=sorted((a.stage58/'cores').glob('*.json'));assert len(fs)==11760
 if a.only_core:fs=[p for p in fs if p.name==a.only_core];assert len(fs)==1
 ih=hashlib.sha256();jobs=[]
 for p in fs:
  raw=p.read_bytes();ih.update(p.name.encode()+b'\0'+raw);jobs.append((p,sha(p)))
 codehash=sha(Path(__file__)); work=[(str(p),q,str(a.out),ih.hexdigest(),codehash,a.max_seconds) for p,q in jobs];res=[]
 with cf.ProcessPoolExecutor(max_workers=a.workers) as ex:
  for z in ex.map(one,work):res.append(z)
 ct=Counter();[ct.update(z['counts']) for z in res]
 atomic(a.out/'summary.json',{'status':'complete','scope':'Necessary root screen of finite stage-5--8 S={2} mothers only; pass is not a factor/cofactor certificate.','code_sha256':codehash,'stage58_input_sha256':ih.hexdigest(),'cores':len(fs),'completed':len(res),'mothers':sum(z['mothers'] for z in res),'counts':dict(ct),'workers':a.workers,'python':sys.version})
 print(json.dumps({'status':'complete','cores':len(fs),'mothers':sum(z['mothers'] for z in res),'counts':dict(ct)},sort_keys=True))
if __name__=='__main__':main()
