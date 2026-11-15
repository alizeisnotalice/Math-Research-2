"""Exhaustive 3x3 verification of the singleton-fiber beta formula."""
import json
from collections import defaultdict
from fractions import Fraction
from itertools import permutations
from exact import grid, peaks

n=3
N=n*n
_,adj=grid(n,n)
groups=defaultdict(list)
for labels in permutations(range(1,N+1)):
    pp=peaks(labels,adj)
    if len(pp)!=2:
        continue
    u,v=sorted(pp,key=lambda z:labels[z],reverse=True)
    H=((1<<N)-1)^(1<<v)
    exposed=[z for z in peaks(labels,adj,H) if z!=u]
    M=max(labels[z] for z in range(N) if adj[v]>>z&1)
    b=max((labels[z] for z in exposed),default=0)
    if M==b:
        continue
    j=N-labels[v]
    a=N-M-1
    horder=tuple(sorted((z for z in range(N) if z!=v),key=lambda z:-labels[z]))
    groups[(v,horder)].append((j,a,M-b,Fraction(M-b,labels[v]-1)))

for values in groups.values():
    js={z[0] for z in values};a=values[0][1];D=values[0][2]
    assert len(values)==a and js==set(range(1,a+1))
    assert all(z[1]==a and z[2]==D and z[3]==Fraction(D,N-z[0]-1) for z in values)
    actual=sum(z[3] for z in values)
    expected=D*sum((Fraction(1,N-j-1) for j in range(1,a+1)),Fraction(0))
    assert actual==expected
print(json.dumps({'n':n,'fibers':len(groups),'beta_positive_outputs':sum(map(len,groups.values())),'verified':True}))
