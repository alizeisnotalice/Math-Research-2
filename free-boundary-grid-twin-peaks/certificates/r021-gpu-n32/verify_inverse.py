#!/usr/bin/env python3
import json,pathlib,sys
p=pathlib.Path(sys.argv[1]);x=json.load(open(p));n=x['n'];N=n*n;P=x['P']
assert x['implementation']=='r021-n32-rank-parallel-terminal-inverse-fibre'
assert x['base_source_sha256']=='6aa1f829dbe97c141781c7440c0f99ef62fa966cbe2c665fab7432a48b2ec9f8'
assert x['method']=='strat' and 2<=n<=32 and P>0
def e(a,b):return abs(a%n-b%n)+abs(a//n-b//n)==1
def d(a,b):return abs(a%n-b%n)+abs(a//n-b//n)
def can(u,v):
 q=[u]
 while q[-1]!=v:
  z=q[-1];q.append(min(w for w in range(N) if e(z,w) and d(w,v)==d(z,v)-1))
 return tuple(q)
def inv(s,dd):
 if d(s[0],s[dd-1])!=dd-1 or any(not e(s[k],s[k+1]) for k in range(dd-1)):return 0
 u=s[0];pp=tuple(s[1:dd]);tau=(u,)+tuple(s[dd:]);ans=0
 for v in range(N):
  if v in s[:dd] or d(u,v)!=dd or not e(pp[-1],v) or can(u,v)!=(u,)+pp+(v,):continue
  dp={(1,0):1}
  for _ in range(N-1):
   nd={}
   for (i,m),z in dp.items():
    placed=set(tau[:i])|{pp[j] for j in range(dd-1) if m>>j&1}
    for j,w in enumerate(pp):
     if not(m>>j&1) and any(e(w,h) for h in placed):nd[i,m|1<<j]=nd.get((i,m|1<<j),0)+z
    if i<len(tau):
     w=tau[i];born=not any(e(w,h) for h in placed)
     if (w==v and born) or (w!=v and not born):nd[i+1,m]=nd.get((i+1,m),0)+z
   dp=nd
  ans+=dp.get((len(tau),(1<<(dd-1))-1),0)
 return ans
r=x['rows'][-1];assert r['terminal'] is True
rep=next(json.loads(line) for line in p.with_suffix('').with_suffix('').read_text().splitlines() if '"kind":"replicate"' in line)
if x['tau']!=N+1:
 assert r['inverse_m2_groups']==['0']*((P+255)//256)
 assert r['inverse_m3_groups']==['0']*((P+255)//256)
 assert all(rep[k] is None for k in ('Zmean_m2','Zmean_m3','mean_m2','mean_m3','scaled_mean_m2','scaled_mean_m3'))
 assert rep['K2'] is None
 print(json.dumps({'kind':'host-terminal-inverse-reference','trace':str(p),'n':n,'P':P,'active':False,'groups_zero':True}))
 raise SystemExit(0)
g2=[];g3=[]
for lo in range(0,P,256):
 a=b=0
 for i in range(lo,min(lo+256,P)):
  w=r['prefix'][i];assert len(w)==N-1 and len(set(w))==N-1 and all(0<=v<N for v in w)
  last=next(v for v in range(N) if v not in set(w));s=tuple(w+[last]);a+=r['allowed'][i]*inv(s,2);b+=r['allowed'][i]*inv(s,3)
 g2.append(a);g3.append(b)
assert r['inverse_m2_groups']==[str(z) for z in g2]
assert r['inverse_m3_groups']==[str(z) for z in g3]
def close(a,b):
 a=float(a);b=float(b)
 return a==b if a==0 or b==0 else abs(a-b)<=5e-12*max(abs(a),abs(b))
for d,groups in ((2,g2),(3,g3)):
 mean=sum(groups)/P
 assert close(rep[f'mean_m{d}'],mean)
 assert close(rep[f'scaled_mean_m{d}'],mean)
 # A finite logZ with Z=0 is a representable non-dead floating underflow.
 # In that case the unscaled product is correctly zero while the scaled
 # terminal mean remains exactly available for a log-domain aggregator.
 assert rep['dead'] is False and rep['logZ'] is not None
 if float(rep['Z'])==0.0:
  assert float(rep[f'Zmean_m{d}'])==0.0
 else:
  assert close(rep[f'Zmean_m{d}'],float(rep['Z'])*mean)
out={'kind':'host-terminal-inverse-reference','trace':str(p),'n':n,'P':P,'m2_groups':[str(z) for z in g2],'m3_groups':[str(z) for z in g3],'m2_total':str(sum(g2)),'m3_total':str(sum(g3)),'overflow':False}
q=p.with_suffix('').with_suffix('.inverse-host.json');q.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
