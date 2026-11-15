"""Independent audit of the proposed beta=0 independent-set recurrence on 3x3."""
from itertools import permutations
import math
import sys
import json
from exact import grid, peaks, vertices

n=int(sys.argv[1]) if len(sys.argv)>1 else 3
assert 2<=n<=4
N=n*n
xy,adj=grid(n,n)
full=(1<<N)-1
closed=[0]*(1<<N)
for I in range(1<<N):
    z=I
    for i in vertices(I):
        z|=adj[i]
    closed[I]=z

# C_G(I): count of labels with every independent I a peak.
C=[0]*(1<<N)
C[0]=math.factorial(N)
for k in range(1,N+1):
    for I in range(1,1<<N):
        if I.bit_count()!=k or any(adj[i]&I for i in vertices(I)):
            continue
        num=sum(C[I^(1<<i)] for i in vertices(I))
        q=closed[I].bit_count()
        assert num%q==0
        C[I]=num//q

def recurrence(v,x):
    T=closed[1<<v]|closed[1<<x]
    U=full^T
    masks=[I for I in range(1<<N) if I&~U==0 and not any(adj[i]&I for i in vertices(I))]
    masks.sort(key=int.bit_count)
    X={};E={}
    for I in masks:
        a=(closed[I]|T).bit_count()
        sx=sum(X[I^(1<<i)] for i in vertices(I))
        se=sum(E[I^(1<<i)] for i in vertices(I))
        assert (C[I]+sx)%(a-1)==0
        X[I]=(C[I]+sx)//(a-1)
        assert X[I]%N==0
        assert (X[I]+se)%a==0
        E[I]=(X[I]+se)//a
    countE=sum((1 if I.bit_count()%2 else -1)*I.bit_count()*E[I] for I in masks if I)
    countVhigh=sum((1 if I.bit_count()%2 else -1)*I.bit_count()*(X[I]//N) for I in masks if I)
    return countE,countVhigh,countE-countVhigh

direct={(v,x):0 for v in range(N) for x in vertices(adj[v])}
if n==3:
    beta0=0
    for labels in permutations(range(1,N+1)):
        pp=peaks(labels,adj)
        if len(pp)!=2:
            continue
        u,v=sorted(pp,key=lambda z:labels[z],reverse=True)
        H=full^(1<<v)
        exposed=[z for z in peaks(labels,adj,H) if z!=u]
        M=max(labels[z] for z in vertices(adj[v]))
        b=max((labels[z] for z in exposed),default=0)
        if b!=M:
            continue
        x=next(z for z in vertices(adj[v]) if labels[z]==M)
        beta0+=1
        T=closed[1<<v]|closed[1<<x]
        vals=sorted((labels[z] for z in vertices(T)),reverse=True)
        assert labels[v]==vals[0] and labels[x]==vals[1]
        direct[(v,x)]+=1
total=0
for edge,expected in direct.items():
    gotE,gotHigh,got=recurrence(*edge)
    if n==3:
        assert got==expected,(edge,gotE,gotHigh,got,expected)
    total+=got
if n==3:
    assert total==beta0
    print(json.dumps({'n':n,'beta_zero_total':beta0,'directed_edge_total':total,'edges_checked':len(direct),'verified':True}))
else:
    print(json.dumps({'n':n,'beta_zero_recurrence_total':total,'directed_edges':len(direct)}))
