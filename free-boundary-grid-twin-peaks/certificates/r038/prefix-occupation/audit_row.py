from fractions import Fraction as F
import json
out=[]
for n in [4,5,8,16,32,64]:
 N=n*n;r=0;nb=lambda v:[Y*n+X for X,Y in [(v%n-1,v//n),(v%n+1,v//n),(v%n,v//n-1),(v%n,v//n+1)] if 0<=X<n and 0<=Y<n]
 pos=list(range(1,N+1));early=[[u for u in nb(v) if pos[u]<pos[v]] for v in range(N)];assert all(early[v] for v in range(1,N))
 gamma=F();psi=F();groups={'top':F(),'column':F(),'interior':F()};counts={k:0 for k in groups}
 for v in range(N):
  if v==r or v in nb(r):continue
  x,y=v%n,v//n;l=min(pos[u] for u in nb(v));b=min((pos[u] for u in nb(v) if early[u]==[v]),default=N+1);R=b-l-1;assert R>0
  if y==0:
   key='top';assert l==x
   if x<=n-2:assert b==x+2 and R==1
  elif x==0:key='column';assert l==(y-1)*n+1 and R==2*n-1
  else:key='interior';assert b==N+1 and R==(n-y+1)*n-x-1 and R>=(n-y)*n
  gamma+=F(1,R);psi+=F(l-1,R);groups[key]+=F(l-1,R);counts[key]+=1
 H=sum((F(1,k) for k in range(1,n)),F());bound=N*(1+H)+1
 assert gamma>=n-3 and psi<=bound and groups['top']<=F(N,2)+1 and groups['column']<=F(N,2) and groups['interior']<=N*H
 out.append(dict(n=n,psi=float(psi),gamma=float(gamma),ratio=float(psi/(N*gamma)),upper_ratio=float((1+H+F(1,N))/(n-3)),counts=counts,all_exact_assertions_pass=True))
print(json.dumps(out,indent=2))
