"""Exact univariate test of the directional area exclusion for rational displacements."""
from itertools import combinations
import json
from sympy import Poly, QQ, Rational as Q, intervals, symbols

s = symbols("s")
den = 1+s*s

def samples(polys):
    # Distinct real roots from all factors are isolated simultaneously.
    ps = [Poly(p,s,domain=QQ) for p in [s, s-1, *polys] if p.expand().has(s)]
    roots = intervals(ps, eps=Q(1,10**20), inf=0, sup=1, strict=True)
    out = [Q(0), Q(1)]
    for ((_,right),_), ((left,_),_) in zip(roots,roots[1:]):
        assert right < left
        out.append((right+left)/2)
    return out

def check(vectors, target_gain, rho=Q(0)):
    assert len(vectors) == 4 and rho >= 0
    T, width, c = Q(".837295006222"), Q("1.245477"), Q(3,1250)
    budget = target_gain + 4*T*rho + Q(22,7)*rho**2
    for sign in (1,-1):
        x,y = sign*(1-s*s),2*s
        along = [a*x+b*y for a,b in vectors]
        across = [0*s]+[-a*y+b*x for a,b in vectors]
        events = along+[a-b for a,b in combinations(across,2)]
        def numerator(t):
            first = sum(p if p.subs(s,t)>=0 else -p for p in along)
            hi = max(across,key=lambda p:p.subs(s,t))
            lo = min(across,key=lambda p:p.subs(s,t))
            return (c*(sign*x+y)*first-(width+2*rho)*(hi-lo)*den-budget*den**2).expand()
        pieces = {numerator(t) for t in samples(events)}
        for t in samples(events+list(pieces)):
            if numerator(t).subs(s,t)>0:
                return dict(sign=sign,parameter=str(t),certified=True)
    return dict(certified=False)

if __name__ == "__main__":
    direction = Q(3,5), Q(4,5)
    moves = [(a*direction[0]-b*direction[1],a*direction[1]+b*direction[0])
             for a,b in zip(map(Q,[".0008","-.0005",".0004","-.0003"]),
                            map(Q,[".00000011","-.00000003",".00000019","0"]))]
    # M*+10^-6 is above the low-area threshold; rho also allows perturbed rotations.
    cases = [check(moves,Q(".000001"),Q(".000000163")),
             check([(Q(0),Q(0))]*4,Q(".000001"))]
    assert [q["certified"] for q in cases] == [True,False]
    print(json.dumps(cases,indent=2))
