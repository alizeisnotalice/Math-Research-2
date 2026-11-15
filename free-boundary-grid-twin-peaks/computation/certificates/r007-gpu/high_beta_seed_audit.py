"""Exact CPU audit for the high-beta, two-peak path-prefix family."""
import json, sys

def nb(v,n):
 x,y=v%n,v//n
 return [z for z,ok in ((v-1,x>0),(v+1,x+1<n),(v-n,y>0),(v+n,y+1<n)) if ok]
def peaks(a,n): return [v for v in range(n*n) if all(a[w]<a[v] for w in nb(v,n))]
def beta(a,n,p):
 v=min(p,key=a.__getitem__);M=b=0
 for x in nb(v,n):
  M=max(M,a[x])
  if all(z==v or a[z]<a[x] for z in nb(x,n)): b=max(b,a[x])
 return (M-b)/(a[v]-1),M,b,a[v]
def seed(n):
 N=n*n;v=N-1;y=(n-2)*n+n-1;order=[0];used={0}
 for r in range(1,n):order.append(r*n);used.add(r*n)
 for c in range(1,n-1):order.append((n-1)*n+c);used.add((n-1)*n+c)
 # Deterministic connected completion avoiding the other v-neighbour y.
 while len(order)<N-2:
  w=next(w for w in range(N) if w not in used and w not in (v,y) and any(z in used for z in nb(w,n)))
  order.append(w);used.add(w)
 order.append(y);a=[0]*N;a[0]=N;a[v]=N-1
 for i,w in enumerate(order[1:],1):a[w]=N-1-i
 return a
for n in map(int,sys.argv[1:] or [8,12,16,24,32]):
 a=seed(n);p=peaks(a,n);d=abs(p[0]%n-p[1]%n)+abs(p[0]//n-p[1]//n);z,M,b,t=beta(a,n,p);N=n*n
 assert len(p)==2 and d==2*(n-1) and (M,b,t)==(N-2*n+2,0,N-1)
 assert z==(N-2*n+2)/(N-2)
 print(json.dumps({'n':n,'K':2,'distance':d,'beta':z,'M':M,'b':b,'t':t,'verified':True}))
