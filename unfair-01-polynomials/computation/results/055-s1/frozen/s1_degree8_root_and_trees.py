#!/usr/bin/env python3
"""Independent root necessary screen and paired prefix trees for degree-8 S={1}.

Input records are only used to identify necessary mother polynomials.  Both
trees restart from G=Q=1, retain signed integer choices and rational 0/1
branches, force only G6=Q6=0 and A7=B7=0, and never treat an unknown tree as
an exclusion.
"""
from __future__ import annotations
import argparse,concurrent.futures as cf,hashlib,json,math,os,time,sys
from collections import Counter
from pathlib import Path

Pair=tuple[int,int]; ZERO=(0,0); ONE=(1,0)
POINTS=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))

def add(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return (-x[0],-x[1])
def scale(k,x):return (k*x[0],k*x[1])
def mul(x,y,T,N):
 p,q=x;r,s=y;return (p*r-q*s*N,p*s+q*r+q*s*T)
def conj(x,T):
 p,q=x;return (p+q*T,-q)
def sign(x,T,D):
 p,q=x;r=2*p+q*T
 if q==0:return (r>0)-(r<0)
 if r>=0 and q>0:return 1
 if r<=0 and q<0:return -1
 left=r*r;right=q*q*D
 return ((left>right)-(left<right)) if r>0 else ((right>left)-(right<left))
def floor_pair(x,T,D):
 p,q=x;B=abs(p)+abs(q)*(abs(T)+math.isqrt(D)+2)+3;lo,hi=-B,B
 while hi-lo>1:
  m=(lo+hi)//2
  if sign((p-m,q),T,D)>=0:lo=m
  else:hi=m
 return lo
def shifts(base,T,D):
 p,q=base
 if q==0:return [(-p,ZERO),(1-p,ONE)]
 h=-floor_pair(base,T,D);z=(p+h,q)
 assert sign(z,T,D)>0 and sign(add(ONE,neg(z)),T,D)>0
 return [(h,z)]
def forced_raw(base,h,T,D):
 z=add(base,(h,0))
 return [(h,z)] if sign(z,T,D)>=0 and sign(add(ONE,neg(z)),T,D)>=0 else []
def forced_zero(base,T,D):
 # The requested A_k=0 is possible exactly when base is an integer.
 if base[1]!=0:return []
 return [(-base[0],ZERO)]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,x):
 t=p.with_suffix(p.suffix+'.partial')
 with t.open('w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)

def conv(X,Y,T,N):
 z=[ZERO]*(len(X)+len(Y)-1)
 for i,u in enumerate(X):
  for j,v in enumerate(Y):z[i+j]=add(z[i+j],mul(u,v,T,N))
 return z
def mother(row):
 T,N,b,c,d,e,f,g=row;D=T*T-4*N;assert D>0 and math.isqrt(D)**2!=D
 P=[ONE,(0,1),(b,0),(c,0),(d,0),(e,0),(f,0),(g,0),ONE]
 Hpair=conv(P,[conj(x,T) for x in P],T,N);assert len(Hpair)==17 and all(q==0 for p,q in Hpair)
 return D,P,[p for p,q in Hpair]
def schur(d):
 """Return first failed degree; None iff all roots strictly in unit disk."""
 while len(d)>1:
  g=math.gcd(*d)
  if d[0]<0:g=-g
  d=[x//g for x in d]
  if abs(d[-1])>=d[0]:return len(d)-1
  first,last=d[0],d[-1];n=len(d)-1
  d=[first*d[i]-last*d[n-i] for i in range(n)]
 return None
def root_screen(row):
 D,P,H=mother(row)
 for j,(p,q) in enumerate(POINTS):
  val=sum(H[i]*p**i*q**(16-i) for i in range(17))
  if val<=0:return {'status':'rejected','code':10+j}
 outer=schur([(v<<i) for i,v in enumerate(H)][::-1])
 if outer is not None:return {'status':'rejected','code':30+outer}
 inner=schur([H[j]<<(16-j) for j in range(17)])
 if inner is not None:return {'status':'rejected','code':50+inner}
 return {'status':'passed','code':255,'H':H}

def state_hash(states):
 h=hashlib.sha256()
 for s in sorted(states,key=lambda z:(z[0],z[1])):h.update(json.dumps([s[0],s[1]],separators=(',',':')).encode()+b'\n')
 return h.hexdigest()
def choices(base,k,T,D,kind):
 if k==6:return forced_raw(base,0,T,D)
 if k==7:return forced_zero(base,T,D)
 return shifts(base,T,D)
def tree_hgq(row,depth,cap,seconds):
 T,N,*_=row;D,P,H=mother(row);states=[((1,),(1,),(1,))];tr=[];begin=time.monotonic();status='unknown_depth'
 for k in range(1,depth+1):
  nxt=[]
  for G,Q,W in states:
   U=ZERO;V=ZERO
   for j in range(1,min(8,k)+1):U=add(U,scale(G[k-j],P[j]));V=add(V,scale(Q[k-j],conj(P[j],T)))
   internal=sum(G[i]*Q[k-i] for i in range(1,k))
   R=internal+sum(H[j]*W[k-j] for j in range(1,min(16,k)+1))
   for g,_Ag in choices(U,k,T,D,'G'):
    for q,_Bq in choices(V,k,T,D,'Q'):
     if R+g+q not in (0,1):continue
     nxt.append((G+(g,),Q+(q,),W+(internal+g+q,)))
  assert len({(x[0],x[1]) for x in nxt})==len(nxt)
  states=nxt;tr.append({'k':k,'live':len(states),'state_sha256':state_hash(states)})
  if not states:status='dead_all_prefixes';break
  if len(states)>cap:status='unknown_state_cap';break
  if time.monotonic()-begin>seconds:status='unknown_time_limit';break
 out={'status':status,'depth_reached':len(tr),'trace':tr}
 if status!='dead_all_prefixes':out['live_prefixes']=[{'G':list(g),'Q':list(q)} for g,q,w in states]
 return out
def coeff_base(P,G,k,T,N):
 z=ZERO
 for j in range(1,min(8,k)+1):z=add(z,scale(G[k-j],P[j]))
 return z
def c_binary(A,B,k,T,N):
 z=ZERO
 for i in range(k+1):z=add(z,mul(A[i],B[k-i],T,N))
 return z[0] if z[1]==0 and z[0] in (0,1) else None
def tree_ab(row,depth,cap,seconds):
 T,N,*_=row;D,P,_=mother(row);Pb=[conj(x,T) for x in P];states=[((1,),(1,),(ONE,),(ONE,))];tr=[];begin=time.monotonic();status='unknown_depth'
 for k in range(1,depth+1):
  nxt=[]
  for G,Q,A,B in states:
   U=coeff_base(P,G,k,T,N);V=coeff_base(Pb,Q,k,T,N)
   for g,Ak in choices(U,k,T,D,'G'):
    for q,Bk in choices(V,k,T,D,'Q'):
     if c_binary(A+(Ak,),B+(Bk,),k,T,N) is None:continue
     nxt.append((G+(g,),Q+(q,),A+(Ak,),B+(Bk,)))
  assert len({(x[0],x[1]) for x in nxt})==len(nxt)
  states=nxt;tr.append({'k':k,'live':len(states),'state_sha256':state_hash(states)})
  if not states:status='dead_all_prefixes';break
  if len(states)>cap:status='unknown_state_cap';break
  if time.monotonic()-begin>seconds:status='unknown_time_limit';break
 out={'status':status,'depth_reached':len(tr),'trace':tr}
 if status!='dead_all_prefixes':out['live_prefixes']=[{'G':list(g),'Q':list(q)} for g,q,A,B in states]
 return out

def one(row,outdir,codehash,sourcehash,depth,cap,seconds):
 target=Path(outdir)/f'mother-{hashlib.sha256(json.dumps(row,separators=(",",":")).encode()).hexdigest()}.json'
 bind={'code_sha256':codehash,'source_sha256':sourcehash,'row':row,'depth':depth,'cap':cap,'seconds_per_tree':seconds}
 if target.exists():
  z=json.loads(target.read_text());assert z['binding']==bind and z['status']=='complete';return z
 screen=root_screen(row)
 if screen['status']!='passed':
  z={'status':'complete','binding':bind,'screen':screen};atomic(target,z);return z
 h=tree_hgq(row,depth,cap,seconds);a=tree_ab(row,depth,cap,seconds)
 assert h['trace']==a['trace']
 z={'status':'complete','binding':bind,'screen':{'status':'passed','code':255,'H':screen['H']},'hgq':h,'ab':a};atomic(target,z);return z

def load_mothers(stage):
 s=json.loads((stage/'summary.json').read_text());assert s['status']=='complete' and s['cores']==s['completed']==1140 and s['row7_states']==67
 fs=sorted((stage/'cores').glob('*.json'));assert len(fs)==1140
 H=hashlib.sha256();rows=set();records=0
 for f in fs:
  raw=f.read_bytes();H.update(f.name.encode()+b'\0'+raw);x=json.loads(raw);assert x['status']=='complete';records+=len(x['records'])
  for b,c,d,e,ff,g,*_ in x['records']:
   q=x['summary'];rows.add((q['T'],q['N'],b,c,d,e,ff,g))
 assert records==67
 return sorted(rows),hashlib.sha256((stage/'summary.json').read_bytes()+H.digest()).hexdigest()

def symbolic_check():
 import sympy as sp
 a,T,N=sp.symbols('a T N');b,c,d,e,f,g=sp.symbols('b c d e f g');G1,G2,G3,G4,G5,G6,G7=sp.symbols('G1:8')
 P=[1,a,b,c,d,e,f,g,1];Pb=[1,T-a,b,c,d,e,f,g,1]
 H=[sp.expand(sum(P[i]*Pb[k-i] for i in range(max(0,k-8),min(8,k)+1))) for k in range(17)]
 assert all(sp.rem(sp.Poly(x,a),sp.Poly(a*a-T*a+N,a)).degree()==0 for x in H)
 A7=sp.expand(sum(P[i]*[1,G1,G2,G3,G4,G5,G6,G7][7-i] for i in range(8)))
 assert sp.expand(A7-(G7+a*G6+b*G5+c*G4+d*G3+e*G2+f*G1+g))==0

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stage57',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=1);ap.add_argument('--depth',type=int,default=128);ap.add_argument('--cap',type=int,default=25000);ap.add_argument('--tree-seconds',type=float,default=600);ap.add_argument('--seconds',type=float,default=14400);ap.add_argument('--mother-limit',type=int);ap.add_argument('--symbolic-check',action='store_true');args=ap.parse_args()
 if args.symbolic_check:symbolic_check();print(json.dumps({'status':'symbolic_check_passed'}));return
 if not __debug__:raise RuntimeError('refuse python -O')
 assert 1<=args.workers<=4
 args.out.mkdir(parents=True,exist_ok=True);rows,sourcehash=load_mothers(args.stage57)
 if args.mother_limit is not None:rows=rows[:args.mother_limit]
 codehash=sha(Path(__file__));start=time.monotonic();done=[];counts=Counter();deadline=start+args.seconds
 with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
  fs=[pool.submit(one,r,str(args.out),codehash,sourcehash,args.depth,args.cap,args.tree_seconds) for r in rows]
  for f in cf.as_completed(fs):
   z=f.result();done.append(z['binding']['row']);counts[z['screen']['status']]+=1
   if z['screen']['status']=='passed':counts['root_pass']+=1;counts['tree_'+z['hgq']['status']]+=1
   print(json.dumps({'done':len(done),'expected':len(rows),'screen':z['screen']['status'],'row':z['binding']['row']}),flush=True)
   if time.monotonic()>deadline:raise TimeoutError('global cap')
 assert len(done)==len(rows)
 out={'status':'complete' if len(rows)==len(load_mothers(args.stage57)[0]) else 'partial_scope','scope':'Necessary degree8 S={1} root screen and prefix trees only; original conjecture remains open.','code_sha256':codehash,'stage57_input_sha256':sourcehash,'input_records':67,'unique_mothers':len(rows),'completed':len(done),'counts':dict(counts),'depth':args.depth,'cap':args.cap,'tree_seconds':args.tree_seconds,'workers':args.workers,'seconds':time.monotonic()-start,'python':sys.version}
 atomic(args.out/'summary.json',out);print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
