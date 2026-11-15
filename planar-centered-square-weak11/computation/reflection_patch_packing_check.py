#!/usr/bin/env python3
"""Fraction certificates for conditional multiple reflection patches (new work after R30)."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,platform

def measure_union(patches):
    edges=sorted({v for p in patches for interval in p['E'] for v in interval})
    area=F(0);q=0;rows=[]
    for left,right in zip(edges,edges[1:]):
        mid=(left+right)/2
        active=[p for p in patches if any(a<mid<b for a,b in p['E'])]
        width=max((p['w'] for p in active),default=F(0))
        area+=2*width*(right-left);q=max(q,len(active))
        if active:rows.append((str(left),str(right),str(width),len(active)))
    raw=sum((2*p['w']*sum((b-a for a,b in p['E']),F(0)) for p in patches),F(0))
    return area,raw,q,rows

def outside(interval,cover):
    pieces=[interval]
    for c,d in cover:
        nxt=[]
        for a,b in pieces:
            if b<=c or d<=a:nxt.append((a,b));continue
            if a<c:nxt.append((a,min(b,c)))
            if d<b:nxt.append((max(a,d),b))
        pieces=nxt
    return pieces

rows=[]
for N in [1,2,5,20,64]:
    # N mass-4 sources at (0,3j), lambda=1. All non-singleton old rectangles empty.
    for n in range(2,N+1):
        assert 9*(n-1)**2>4*n
    w=F(1,5);height=F(2,5);rho=F(7,5)
    assert rho*rho<2 and F(6,5)<rho
    avg=F(8)/(4*rho*rho);assert avg==F(50,49)>1
    old_vertical=[(F(3*j-1),F(3*j+1)) for j in range(N)]
    patches=[]
    for j in range(N):
        y=F(3*j)
        E=[(y-F(7,5),y-1),(y+1,y+F(7,5))]
        assert all(outside(interval,old_vertical)==[interval] for interval in E)
        patches.append(dict(w=w,E=E,source=[j]))
    area,raw,q,_=measure_union(patches)
    assert area==raw==F(8*N,25) and q==1
    M=F(4*N);gain=area/(2*M);assert gain==F(1,25)
    D=16*N+sum((2*(2*p['w']+(b-a)) for p in patches for a,b in p['E']),F(0))
    assert D==F(96*N,5)
    eta=area/(2*D);assert eta==F(1,120)
    smooth_gain=(area-D*eta)/(2*M);assert smooth_gain==F(1,50)
    rows.append(dict(N=N,area=str(area),receiver_multiplicity=q,covered_source_mass=str(M),gain=str(gain),smoothing_eta=str(eta),smooth_gain=str(smooth_gain),mixed_average=str(avg)))
# Exact centered-strip overlap bookkeeping: area 5, raw sum 6, overlap cost 1.
A,raw,q,cells=measure_union([dict(w=F(1),E=[(F(0),F(2))]),dict(w=F(1,2),E=[(F(1),F(3))])])
assert A==5 and raw==6 and raw-A==1 and q==2 and A>=raw/q
# Small source at (0,0), mass4; large source at (-1,0), mass12.
# Pair radius2 has right edge1, vertical interval(-2,2); small raw strips are hidden.
assert F(3)<4 and F(2)<4  # sqrt3<2 and sqrt2<2, with positive compared quantities
hidden=[(F(-7,5),F(-1)),(F(1),F(7,5))]
assert all(outside(E,[(F(-2),F(2))])==[] for E in hidden)
out=dict(python=platform.python_version(),arithmetic='Fraction exact rational certificates',seed=None,tolerance=None,examples=rows,overlap=dict(union_area=str(A),raw_area=str(raw),cost=str(raw-A),q=q,cells=cells),shielded_small_witness='PASS: rational top/bottom patches fully covered by pair old interval',status='PASS',scope='finite certificates only; all-N and general patch theorem proved in separate report',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
print(json.dumps(out,indent=2))
