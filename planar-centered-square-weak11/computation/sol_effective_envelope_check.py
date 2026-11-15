from fractions import Fraction as F
from itertools import combinations,permutations
from pathlib import Path
import json,platform,hashlib,math

def area(rects):
    xs=sorted({x for z in rects for x in z[:2]});total=F(0)
    for l,r in zip(xs,xs[1:]):
        intervals=sorted((b,t) for a,c,b,t in rects if a<r and c>l)
        covered=F(0);end=None
        for b,t in intervals:
            if end is None or b>end:covered+=t-b;end=t
            elif t>end:covered+=t-end;end=t
        total+=(r-l)*covered
    return total

a=F(49,100);points=[(x,y) for x in (-a,a) for y in (-a,a)];mass=F(1,4)
def span(indices):
    q=[points[i] for i in indices]
    return max(max(z[k] for z in q)-min(z[k] for z in q) for k in (0,1))
def square(i,r):
    x,y=points[i];return (x-r,x+r,y-r,y+r)

subsets=[];realrects=[]
for k in range(1,5):
    for S in combinations(range(4),k):
        m=k*mass;d=span(S);valid=d*d<m
        subsets.append(dict(indices=list(S),mass=str(m),span=str(d),span_squared=str(d*d),eligible=valid))
        assert valid==(k in (1,4))
        if k==1:realrects.append(square(S[0],F(1,4)))
        if k==4:realrects.append((-F(1,100),F(1,100),-F(1,100),F(1,100)))
assert area(realrects)==F(2501,2500)
orders=[]
for order in permutations(range(4)):
    old=[];eff=[];rects=[]
    for i in order:
        value=F(0)
        for k in range(len(old)+1):
            for S in combinations(old,k):
                m=k*mass;d=span(S+(i,))
                if d*d<m+mass:value=max(value,m)
        eff.append(value)
        assert value in (F(0),F(3,4))
        r=F(1,4) if value==0 else F(1,2)
        rects.append(square(i,r));old.append(i)
    assert eff==[F(0),F(0),F(0),F(3,4)]
    u=area(rects);assert u==F(7,4)
    orders.append(dict(order=list(order),effective_masses=[str(x) for x in eff],union_area=str(u)))
proof=Path(__file__).resolve().parents[1]/'proofs/sol-effective-envelope.md'
sixpoints=[(x,y) for x in (-a,F(0),a) for y in (-a,a)]
sixchecks=[]
for k in range(1,7):
    for S in combinations(range(6),k):
        pts=[sixpoints[i] for i in S]
        d=max(max(z[q] for z in pts)-min(z[q] for z in pts) for q in (0,1))
        eligible=d*d<F(k,6)
        same_row=len({z[1] for z in pts})==1
        expected=k==1 or k==6 or (k==2 and same_row and d==a)
        assert eligible==expected
        sixchecks.append(dict(indices=list(S),mass=str(F(k,6)),span_squared=str(d*d),eligible=eligible))
assert F(3,2)**2>2 and F(2,5)**2<F(1,6)
assert F(5,3)-F(1,4)+F(48,125)-F(7,4)==F(19,375)
print(json.dumps(dict(python=platform.python_version(),method='exact Fraction full four-point subset/order enumeration and six-point subset classification; exact rectangular union; single process, no random seed or tolerance',a=str(a),masses=['1/4']*4,subsets=subsets,orders=orders,true_halo_area='2501/2500',min_envelope_area='7/4',six_point_subsets=sixchecks,six_point_minimum_formula='5/3-sqrt(2)/6+(24/25)/sqrt(6)',six_point_float_display=5/3-math.sqrt(2)/6+(24/25)/math.sqrt(6),six_point_above_7_over_4_by_more_than='19/375',proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest()),indent=2))
