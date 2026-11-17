#!/usr/bin/env python3
import sys, math
from fractions import Fraction
p=open(sys.argv[1]).read().splitlines(); assert len(p)==2 and p[1].startswith('word ')
h=dict(zip(p[0].split()[::2],p[0].split()[1::2])); w=list(map(int,p[1].split()[1:])); n=int(h['n']);N=n*n;r=int(h['root'])
assert n>=2 and len(w)==N and 0<=r<N and set(w)==set(range(N)) and w[0]==r
pos=[0]*N
for i,x in enumerate(w):pos[x]=i+1
nb=lambda v:[Y*n+X for X,Y in ((v%n-1,v//n),(v%n+1,v//n),(v%n,v//n-1),(v%n,v//n+1)) if 0<=X<n and 0<=Y<n]
assert w[0]==r and all(any(pos[u]<pos[v] for u in nb(v)) for v in range(N) if v!=r)
ans=Fraction()
for x in range(N):
 if x==r or x in nb(r):continue
 l=min(pos[u] for u in nb(x)); b=min((pos[v] for v in nb(x) if [u for u in nb(v) if pos[u]<pos[v]]==[x]),default=N+1); R=b-l-1; assert R>0; ans+=Fraction(l-1,R)
got=float(h['psi']); norm=float(h['norm']); exact=float(ans); target=exact/(N*math.log(n)); tol=1e-11*max(1.0,abs(exact))
assert math.isfinite(got) and math.isfinite(norm) and abs(got-exact)<=tol and abs(norm-target)<=1e-11*max(1.0,abs(target))
print(ans, exact)
