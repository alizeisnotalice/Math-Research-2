#!/usr/bin/env python3
import sys
from pathlib import Path
mode,n,path=sys.argv[1],int(sys.argv[2]),Path(sys.argv[3]);N=n*n;A=[set() for _ in range(N)]
for y in range(n):
 for x in range(n):
  k=y*n+x
  for X,Y in ((x+1,y),(x-1,y),(x,y+1),(x,y-1)):
   if 0<=X<n and 0<=Y<n:A[k].add(Y*n+X)
def ok(w,t):
 s=set();b=[]
 for k,v in enumerate(w,1):
  if not(A[v]&s):b.append(k)
  s.add(v)
 return b==([1] if len(w)<t else [1,t])
c=0
for q in path.read_text().splitlines():
 z=q.split(','); pre=list(map(int,z[1].split(':')));t,r=int(z[2]),int(z[3]);assert len(pre)==r and len(set(pre))==r and all(0<=v<N for v in pre) and 2<=t<=N and ok(pre,t)
 if mode=='b1':
  i,v,accept=int(z[9]),int(z[10]),int(z[11]);post=list(map(int,z[12].split(':')))
  assert 0<=i<r and 0<=v<N
  ri,xi,rv,xv=map(int,z[5:9]);assert xi>=0 and xv>=0 and ri>=((1<<64)%r) and rv>=((1<<64)%N) and ri%r==i and rv%N==v
  y=pre[:];j=y.index(v) if v in y else -1
  if j>=0:y[i],y[j]=y[j],y[i]
  else:y[i]=v
  assert accept==ok(y,t) and post==(y if accept else pre)
 else:
  i1,v1,i2,v2,accept=map(int,z[5:10]);post=list(map(int,z[10].split(':')));y=pre[:]
  for i,v in ((i1,v1),(i2,v2)):
   assert 0<=i<r and 0<=v<N
   j=y.index(v) if v in y else -1
   if j>=0:y[i],y[j]=y[j],y[i]
   else:y[i]=v
  assert accept==ok(y,t) and post==(y if accept else pre)
 assert len(post)==r and len(set(post))==r and ok(post,t) and r<N;c+=1
print({'status':'passed','mode':mode,'events':c})
