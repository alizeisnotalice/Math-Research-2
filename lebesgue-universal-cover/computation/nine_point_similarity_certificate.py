"""Rational outer host and four-variable similarity duals for the prescribed nine points."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse
import json
from math import ceil, pi, tan
from flint import arb, ctx
from certify_b9_padding_four_minimizer import Hull, dot, features, q

DEN = 10**12+10
POINTS = [(F(x,DEN),F(s*y,DEN)) for x,y in
          [(-477094447069,173648177667),(-253856652971,439692620786),
           (88163490354,500000000000),(388930956715,326351822333)] for s in (-1,1)]
POINTS.append((F(507713305943,DEN),F(0)))
HERE = Path(__file__).resolve().parent
CERT = HERE/"nine-point-similarity-duals.json"

def directions(count):
    for quadrant in range(4):
        for j in range(count):
            t=F(j,count)
            x,y=(1-t*t)/(1+t*t),2*t/(1+t*t)
            for _ in range(quadrant):
                x,y=-y,x
            yield x,y

def support(host,u):
    # With exact face ties, vertices give valid lower witnesses on either branch.
    v=tuple(map(q,u))
    angle=arb.atan2(v[1],v[0])
    low=high=arb(1)/2
    for n in (3,5,7):
        for k,rad,interior in features(n,angle,host.rot[n]):
            base=dot(v,host.vertices[n][k])
            high=max(high,(base+rad).upper())
            low=max(low,(base+(rad if interior else 0)).lower())
    return low,high

def clip(poly,n,h):
    out=[]
    for a,b in zip(poly,poly[1:]+poly[:1]):
        l,r=dot(a,n)-h,dot(b,n)-h
        if l<=0:
            out.append(a)
        if l*r<0:
            out.append(tuple((r*x-l*y)/(r-l) for x,y in zip(a,b)))
    return out

def host_polygon():
    host=Hull()
    edges=[]
    normals=list(directions(16))
    for _,angle,*_ in host.pieces:
        # A rational unit normal close to the meeting direction of two pieces.
        quadrant, offset = divmod(float(angle.mid()) % (2*pi),pi/2)
        t=F(round(tan(offset/2)*10**7),10**7)
        x,y=(1-t*t)/(1+t*t),2*t/(1+t*t)
        for _ in range(int(quadrant)):
            x,y=-y,x
        if (x,y) not in normals:
            normals.append((x,y))
    for u in normals:
        _,upper=support(host,u)
        top=F(ceil(float(upper)*10**10)+2,10**10)
        assert q(top)>upper
        edges.append((u,top))
    poly=[(F(x),F(y)) for x,y in ((-2,-2),(2,-2),(2,2),(-2,2))]
    for u,h in edges:
        poly=clip(poly,u,h)
    assert min(len(poly),len(set(poly)))==len(poly) and len(poly)>=3
    assert all(dot(v,v)<1 for v in poly)
    for a,b,c in zip(poly,poly[1:]+poly[:1],poly[2:]+poly[:2]):
        assert (b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0])>0
    # Bound h_Q-h_K on a second rational unit grid, then extend by Lipschitz continuity.
    density=2048
    worst=arb(0)
    k=max(range(len(poly)),key=lambda i:poly[i][0])
    for u in directions(density):
        while dot(poly[(k+1)%len(poly)],u)>dot(poly[k],u):
            k=(k+1)%len(poly)
        assert dot(poly[k-1],u)<=dot(poly[k],u) and dot(poly[(k+1)%len(poly)],u)<=dot(poly[k],u)
        low,_=support(host,u)
        worst=max(worst,(q(dot(poly[k],u))-low).upper())
    error=worst+q(F(2,density))
    assert error<q(F('.00392897')),str(error)
    return edges,poly,error

def rows_for(edges,points):
    return [((x*u+y*v,x*v-y*u,u,v),h) for (u,v),h in edges for x,y in points]

def generate(rows,targets):
    import numpy as np
    from scipy.optimize import linprog
    M=np.array([list(map(float,p)) for p,_ in rows])
    rhs=np.array([float(h) for _,h in rows])
    records=[]
    for u,v in targets:
        ans=linprog([-float(u),-float(v),0,0],A_ub=M,b_ub=rhs,bounds=[(None,None)]*4,method="highs-ds")
        assert ans.success
        records.append([[i,round(-a*10**12)] for i,a in enumerate(ans.ineqlin.marginals) if -a>1e-13])
    return records

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--generate',action='store_true')
    args=p.parse_args()
    ctx.prec=160
    diameter=max(sum((x-y)**2 for x,y in zip(a,b)) for a,b in combinations(POINTS,2))
    assert diameter==F(10**24,DEN**2) and F(9,10)<diameter<1
    assert all(dot(p,p)<1 for p in POINTS)
    edges,poly,error=host_polygon()
    density=12
    targets=list(directions(density))
    data={"weights":[],"host":[[str(x),str(y),str(h)] for (x,y),h in edges]}
    if not args.generate:
        frozen=json.loads(CERT.read_bytes())
        assert frozen["host"]==data["host"] and len(frozen["weights"])==2
    # Positive control: the rational unit square universally covers every diameter-one input.
    # Recover a unit placement (translation aligns both lower extrema) for each hand.
    positive=[]
    for hand in (1,-1):
        pts=[(x,hand*y) for x,y in POINTS]
        tx=-min(x for x,_ in pts)
        ty=-min(y for _,y in pts)
        assert all(0<=x+tx<=1 and 0<=y+ty<=1 for x,y in pts)
        positive.append([str(tx),str(ty)])
    bounds=[]
    terms=0
    for branch in (0,1):
        points=[(x,-y if branch else y) for x,y in POINTS]
        rows=rows_for(edges,points)
        weights=generate(rows,targets) if args.generate else frozen["weights"][branch]
        assert len(weights)==len(targets)
        data["weights"].append(weights)
        top=F(0)
        for terms_for_ray,u in zip(weights,targets):
            target=(*u,F(0),F(0)); coeff=[F(0)]*4; h=F(0)
            for i,v in terms_for_ray:
                assert type(i) is int and 0<=i<len(rows) and type(v) is int and v>0
                w=F(v,10**12);row,b=rows[i]
                for j in range(4):
                    coeff[j]+=w*row[j]
                h+=w*b
            # Q lies inside the unit disk. Its diameter is <=2, and diam(P)>0.9;
            # the similarity scale is <3, so each translation coordinate has modulus <4.
            h+=4*sum(abs(a-b) for a,b in zip(coeff,target))
            top=max(top,h);terms+=len(terms_for_ray)
        bound=top/(1-F(1,2*density*density))
        assert bound<1,str(float(bound))
        bounds.append(str(bound))
    if args.generate:
        CERT.write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(json.dumps(dict(halfplanes=len(edges),polygon_edges=len(poly),distance_upper=str(error),
              directions_per_branch=len(targets),dual_nonzero_terms=terms,
              scale_upper=bounds,scale_upper_display=list(map(lambda t:float(F(t)),bounds)),
              unit_square_translations=positive,diameter_squared=str(diameter),certificate_sha256=sha256(CERT.read_bytes()).hexdigest()),indent=2))

if __name__=='__main__':
    main()
