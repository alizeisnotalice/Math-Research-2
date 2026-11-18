from fractions import Fraction
from pathlib import Path
import json,sys
from audit import coefficients,peaks

def neighbors(n):
    result=[]
    for v in range(n*n):
        x,y=v%n,v//n
        result.append(tuple(Y*n+X for X,Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1))
                            if 0<=X<n and 0<=Y<n))
    return result

def run(n,rootpair,v,nbr):
    # Deterministic descending word with the prescribed nonadjacent peaks.
    N=n*n
    order=list(rootpair)
    used=set(order)
    for u in order:
        for w in nbr[u]:
            if w not in used:
                used.add(w)
                order.append(w)
    assert len(order)==N
    label=[0]*N
    for i,u in enumerate(order):label[u]=N-i
    assert peaks(label,nbr)==tuple(sorted(rootpair))
    assert set(rootpair)<=set(nbr[v])
    c=(v%n+2*(v//n))%5
    l,r=(n-n//4,n) if 2*(v%n)<n else (0,n//4)
    sites=[w for w in range(N) if l<=w%n<r and (w%n+2*(w//n))%5==c]
    blocks=[set((w,*nbr[w])) for w in [v,*sites]]
    occupied=set.union(*blocks)
    assert len(occupied)==sum(map(len,blocks))
    assert all(any(abs(p%n-q%n)+abs(p//n-q//n)>n/8 for q in block) for p in blocks[0] for block in blocks[1:])
    assert all(abs(p%n-q%n)+abs(p//n-q//n)>n/8 for p in blocks[0] for block in blocks[1:] for q in block)
    central=coefficients(label,v,nbr,N+1)
    assert central[0]==0 and central[1]>0 and central[2]>0
    assert central[2]<=N*central[1]
    others=[coefficients(label,w,nbr,N+1) for w in sites]
    assert all(row[0]>0 for row in others)
    high=sum(200*x>199*(N+1) for x in label)
    assert high*100<=N and N*sorted(label)[-2]<=(N-1)*(N+1)
    admissible=sum(200*max(label[w] for w in nbr[z])<=199*(N+1) for z in sites)
    assert admissible*40>=N
    for w,row in zip(sites,others):
        if 200*max(label[z] for z in nbr[w])<=199*(N+1):
            assert 200*row[1]>=row[0]
    mass_floor=Fraction(admissible,200)
    central_odds=Fraction(central[2],central[1])
    assert mass_floor>=central_odds/8000
    return dict(n=n,rootpair=rootpair,common_center=v,far_centers=len(sites),
                admissible=admissible,local_p2_over_p1=str(central_odds))

rows=[]
for n in (128,192,256):
    adj=neighbors(n)
    for p,q,v in (((0,0),(2,0),(1,0)),((n-1,n-2),(n-2,n-1),(n-1,n-1)),
                  ((n//2,n//2),(n//2+1,n//2+1),(n//2,n//2+1))):
        idx=lambda p:p[1]*n+p[0]
        rows.append(run(n,(idx(p),idx(q)),idx(v),adj))
Path(sys.argv[1]).write_text(json.dumps(rows,indent=2)+'\n')
print('9 near-source geometric and integer odds checks')
