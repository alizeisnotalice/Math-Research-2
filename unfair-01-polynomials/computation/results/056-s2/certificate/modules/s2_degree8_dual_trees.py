#!/usr/bin/env python3
"""Independent AB and H(GQ) prefix trees for root-screen S={2} mothers.

Both trees restart at G=Q=1 from the extracted root-screen input.  They impose
only the established G4=Q4=0 and A6=B6=0 constraints, retain every 0/1
rational endpoint, and mark depth/time/state-cap survivors as unknown.
"""
from __future__ import annotations
import argparse,concurrent.futures as cf,gzip,hashlib,json,math,os,sys,time
from collections import Counter
from pathlib import Path
Pair=tuple[int,int]; ZERO=(0,0); ONE=(1,0)
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return (-x[0],-x[1])
def scale(k,x):return (k*x[0],k*x[1])
def mul(x,y,T,N):
 p,q=x;r,s=y;return (p*r-q*s*N,p*s+q*r+q*s*T)
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
def forced_raw(base,h,T,D):
 z=add(base,(h,0))
 return [(h,z)] if sign_plus(z,T,D)>=0 and sign_plus(add(ONE,neg(z)),T,D)>=0 else []
def forced_zero(base,T,D):
 return [(-base[0],ZERO)] if base[1]==0 else []
def choices(base,k,T,D):
 if k==4:return forced_raw(base,0,T,D) # G4=Q4=0, still require A4/B4 in [0,1]
 if k==6:return forced_zero(base,T,D)  # A6=B6=0
 return shifts(base,T,D)
def binary(x):return x[1]==0 and x[0] in (0,1)
def state_hash(states):
 h=hashlib.sha256()
 for z in sorted(states,key=lambda x:(x[0],x[1])):h.update(json.dumps([z[0],z[1]],separators=(',',':')).encode()+b'\n')
 return h.hexdigest()
def mother(row):
 L,T,N,D,b,c,d,e,f=row
 P=[ONE,(L,0),(0,1),(b,0),(c,0),(d,0),(e,0),(f,0),ONE]
 F=[1,L,0,b,c,d,e,f,1];H=[0]*17
 for i,x in enumerate(F):
  for j,y in enumerate(F):H[i+j]+=x*y
 for i,x in enumerate(F):H[i+2]+=T*x
 H[4]+=N
 return P,H
def coeff(P,G,k,T,N):
 z=ZERO
 for j in range(1,min(8,k)+1):z=add(z,scale(G[k-j],P[j]))
 return z
def tree_ab(row,depth,cap,seconds):
 L,T,N,D,*_=row;P,H=mother(row);Pb=[conj(x,T) for x in P];states=[((1,),(1,),(ONE,),(ONE,))];trace=[];begin=time.monotonic();status='unknown_depth'
 for k in range(1,depth+1):
  nxt=[]; rejected=0
  for G,Q,A,B in states:
   U=coeff(P,G,k,T,N);V=coeff(Pb,Q,k,T,N)
   for g,Ak in choices(U,k,T,D):
    for q,Bk in choices(V,k,T,D):
     C=ZERO
     for i in range(k+1):C=add(C,mul((A+(Ak,))[i],(B+(Bk,))[k-i],T,N))
     if not binary(C):rejected+=1;continue
     nxt.append((G+(g,),Q+(q,),A+(Ak,),B+(Bk,)))
  assert len({(z[0],z[1]) for z in nxt})==len(nxt)
  states=nxt;trace.append({'k':k,'live':len(states),'state_sha256':state_hash(states),'nonbinary_C':rejected})
  if not states:status='dead_all_prefixes';break
  if len(states)>cap:status='unknown_state_cap';break
  if time.monotonic()-begin>seconds:status='unknown_time_limit';break
 out={'status':status,'depth_reached':len(trace),'trace':trace}
 if status!='dead_all_prefixes':out['live_prefix_count']=len(states)
 return out
def tree_hgq(row,depth,cap,seconds):
 L,T,N,D,*_=row;P,H=mother(row);states=[((1,),(1,),(1,))];trace=[];begin=time.monotonic();status='unknown_depth'
 for k in range(1,depth+1):
  nxt=[];rejected=0
  for G,Q,W in states:
   U=coeff(P,G,k,T,N);V=coeff([conj(x,T) for x in P],Q,k,T,N)
   internal=sum(G[i]*Q[k-i] for i in range(1,k)); R=internal+sum(H[j]*W[k-j] for j in range(1,min(16,k)+1))
   for g,_ in choices(U,k,T,D):
    for q,_ in choices(V,k,T,D):
     if R+g+q not in (0,1):rejected+=1;continue
     nxt.append((G+(g,),Q+(q,),W+(internal+g+q,)))
  assert len({(z[0],z[1]) for z in nxt})==len(nxt)
  states=nxt;trace.append({'k':k,'live':len(states),'state_sha256':state_hash(states),'nonbinary_C':rejected})
  if not states:status='dead_all_prefixes';break
  if len(states)>cap:status='unknown_state_cap';break
  if time.monotonic()-begin>seconds:status='unknown_time_limit';break
 out={'status':status,'depth_reached':len(trace),'trace':trace}
 if status!='dead_all_prefixes':out['live_prefix_count']=len(states)
 return out
def atomic(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.partial')
 with t.open('w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def one(job):
 row,outdir,inputhash,codehash,depth,cap,seconds=job;rid=hashlib.sha256(json.dumps(row,separators=(',',':')).encode()).hexdigest();target=Path(outdir)/f'mother-{rid}.json'
 bind={'version':1,'row':row,'input_payload_sha256':inputhash,'code_sha256':codehash,'depth':depth,'cap':cap,'seconds_per_tree':seconds}
 if target.exists():
  x=json.loads(target.read_text())
  if x.get('binding')==bind and x.get('status')=='complete':return x
  raise RuntimeError('stale tree cache '+str(target))
 a=tree_ab(row,depth,cap,seconds);h=tree_hgq(row,depth,cap,seconds)
 if a['trace']!=h['trace'] or a['status']!=h['status']:raise AssertionError('AB/HGQ mismatch '+str(row))
 x={'status':'complete','binding':bind,'ab':a,'hgq':h};atomic(target,x);return x
def load_input(p:Path):
 meta=json.loads(p.with_suffix('').with_suffix('.json').read_text()) if p.name.endswith('.jsonl.gz') else None
 with gzip.open(p,'rt') as f:rows=[json.loads(x) for x in f if x.strip()]
 payload=b''.join(json.dumps(x,separators=(',',':')).encode()+b'\n' for x in rows)
 assert len(rows)==len({tuple(x) for x in rows})==meta['mothers'] and rows==sorted(rows)
 assert meta and meta['rows_payload_sha256']==hashlib.sha256(payload).hexdigest()
 return rows,meta['rows_payload_sha256']
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,default=HERE/'s2-tree-input-v1.jsonl.gz');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=1);ap.add_argument('--depth',type=int,default=128);ap.add_argument('--cap',type=int,default=25000);ap.add_argument('--tree-seconds',type=float,default=10.0);ap.add_argument('--limit',type=int);a=ap.parse_args()
 # The local bounded runs use at most four workers; the user-authorized
 # MacBook production batch has thirteen physical CPU workers available.
 assert 1<=a.workers<=13;rows,inputhash=load_input(a.input);full=True
 if a.limit is not None:rows=rows[:a.limit];full=False
 codehash=sha(Path(__file__));start=time.monotonic();done=[];ct=Counter()
 with cf.ProcessPoolExecutor(max_workers=a.workers) as pool:
  for x in pool.map(one,[(r,str(a.out),inputhash,codehash,a.depth,a.cap,a.tree_seconds) for r in rows]):
   done.append(x['binding']['row']);ct['tree_'+x['ab']['status']]+=1
 assert len(done)==len(rows)
 summary={'status':'complete' if full else 'partial_scope','scope':'Independent necessary AB/HGQ prefix trees from the extracted root-screen pass mothers; unknown status is not an exclusion.','code_sha256':codehash,'input_payload_sha256':inputhash,'mothers':len(rows),'completed':len(done),'counts':dict(ct),'depth':a.depth,'cap':a.cap,'tree_seconds':a.tree_seconds,'workers':a.workers,'seconds':time.monotonic()-start,'python':sys.version}
 atomic(a.out/'summary.json',summary);print(json.dumps(summary,sort_keys=True))
if __name__=='__main__':main()
