#!/usr/bin/env python3
import json,pathlib,sys
p=pathlib.Path(sys.argv[1]);x=json.load(open(p));n=x['n'];N=n*n;P=x['P'];assert x['tau']==N+1
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
r=x['rows'][-1];g2=[];g3=[]
for lo in range(0,P,256):
 a=b=0
 for i in range(lo,min(lo+256,P)):
  w=r['prefix'][i];last=next(v for v in range(N) if v not in set(w));s=tuple(w+[last]);a+=r['allowed'][i]*inv(s,2);b+=r['allowed'][i]*inv(s,3)
 g2.append(a);g3.append(b)
out={'kind':'host-terminal-inverse-reference','trace':str(p),'n':n,'P':P,'m2_groups':[str(z) for z in g2],'m3_groups':[str(z) for z in g3],'m2_total':str(sum(g2)),'m3_total':str(sum(g3)),'overflow':False}
q=p.with_suffix('').with_suffix('.inverse-host.json');q.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
