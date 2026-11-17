"""Exact finite implementation checks; general theorems are in proofs/001-bridge.md."""
from fractions import Fraction as F
from itertools import combinations, product
from math import isqrt
import random, json, sys, hashlib

def cross(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def hull(points):
    p=sorted(set(points))
    if len(p)<2:return p
    lo=[];hi=[]
    for q in p:
        while len(lo)>1 and cross(lo[-2],lo[-1],q)<=0:lo.pop()
        lo.append(q)
    for q in reversed(p):
        while len(hi)>1 and cross(hi[-2],hi[-1],q)<=0:hi.pop()
        hi.append(q)
    return lo[:-1]+hi[:-1]
def sqdist(p,q):return sum((a-b)**2 for a,b in zip(p,q))
def exact_sqrt(q):
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    assert a*a==q.numerator and b*b==q.denominator
    return F(a,b)
N=[(F(1),F(0)),(F(3,5),F(4,5)),(F(0),F(1)),(F(-3,5),F(4,5)),(F(-1),F(0)),(F(-3,5),F(-4,5)),(F(0),F(-1)),(F(3,5),F(-4,5))]
rng=random.Random(20261117); checks=[]
for k in range(101):
    offsets=[F(1,2)+(F(rng.randrange(6),10) if k else 0) for _ in N]
    pts=[]
    for i,j in combinations(range(8),2):
        a,b=N[i];c,d=N[j];det=a*d-b*c
        if not det:continue
        x,y=(offsets[i]*d-b*offsets[j])/det,(a*offsets[j]-offsets[i]*c)/det
        if all(u*x+v*y<=w for (u,v),w in zip(N,offsets)):pts.append((x,y))
    H=hull(pts);E=list(zip(H,H[1:]+H[:1]))
    A=sum(p[0]*q[1]-p[1]*q[0] for p,q in E)/2
    P=sum(exact_sqrt(sqdist(p,q)) for p,q in E)
    assert A>=P/4
    assert all(x*x+y*y<=4*A*A for x,y in H)
    if k==0:assert A==P/4
    checks.append({'area':str(A),'perimeter':str(P),'slack':str(A-P/4)})
netcounts={}
for n in (1,2):
    grid=[(F(i,n),F(j,n)) for i,j in product(range(-n,n+1),repeat=2) if i or j]
    count=0
    sizes=range(len(grid)+1) if n==1 else range(3)
    for size in sizes:
        for tail in combinations(grid,size):
            Q=[(F(0),F(0)),*tail]
            if any(sqdist(a,b)>(1+F(2,n))**2 for a,b in combinations(Q,2)):continue
            H=hull([(x/(1+F(2,n)),y/(1+F(2,n))) for x,y in Q])
            assert all(sqdist(a,b)<=1 for a,b in combinations(H,2))
            count+=1
    netcounts[str(n)]={'accepted':count,'scope':'all subsets containing origin' if n==1 else 'at most two nonzero grid points'}
# delta=3/n; pi<4 yields entirely rational certified error envelopes with U=1.
budgets=[]
for eps in (F(1,10),F(1,100),F(1,1000),F(1,10000)):
    n=(24*eps.denominator+eps.numerator-1)//eps.numerator
    delta=F(3,n);error=4*delta+4*delta*delta
    assert error<=eps
    budgets.append({'epsilon':str(eps),'n':n,'delta':str(delta),'error_upper_pi_replaced_by_4':str(error)})
print(json.dumps({'python':sys.version,'script_sha256':hashlib.sha256(open(__file__,'rb').read()).hexdigest(),'seed':20261117,'polygon_checks':len(checks),'tangential_example':checks[0],'finite_grid_checks':netcounts,'error_budgets':budgets,'status':'passed','scope':'finite exact checks, not solution of any finite-family minimum'},indent=2))
