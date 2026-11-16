from itertools import permutations
from fractions import Fraction
import json
from pathlib import Path
n=3;N=9
nb=[set() for _ in range(N)]
for v in range(N):
 x,y=v%n,v//n
 for X,Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)):
  if 0<=X<n and 0<=Y<n:nb[v].add(Y*n+X)
def peaks(w):
 s=set();p=[]
 for x in w:
  if not(nb[x]&s):p.append(x)
  s.add(x)
 return p
cov=[]
for a in range(N):
 for b in nb[a]:
  if a>b:continue
  ca=nb[a]-{b};cb=nb[b]-{a}
  for P in __import__('itertools').combinations(range(N),3):
   S=set(P)
   if any(nb[x]&S-{x} for x in S) or set.intersection(*(nb[x] for x in S)):continue
   if set(P)<=ca|cb and len(S&ca)==2 and len(S&cb)==1:
    uv=sorted(S&ca);q=next(iter(S&cb));cov.append((frozenset(S),a,b,uv[0],uv[1],q))
   if set(P)<=ca|cb and len(S&cb)==2 and len(S&ca)==1:
    uv=sorted(S&cb);q=next(iter(S&ca));cov.append((frozenset(S),b,a,uv[0],uv[1],q))
cnt=[0]*len(cov);gain=[Fraction() for _ in cov]
for w in permutations(range(N)):
 P=frozenset(peaks(w))
 if len(P)!=3:continue
 pos={x:i for i,x in enumerate(w)}
 for t,(Q,a,b,u,v,q) in enumerate(cov):
  if P==Q:
   cnt[t]+=1;k=max(pos[u],pos[v]);gain[t]+=Fraction(abs(pos[u]-pos[v]),N-k-1)
out={"covers":len(cov),"count":sum(cnt),"gain":str(sum(gain,Fraction())),"per":[[cnt[i],str(gain[i])] for i in range(len(cov))]}
assert out["count"]==129024 and out["gain"]=="51072"
Path(__file__).with_name("n3-direct-audit.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
