"""Verify area inequalities for a polygon and an arbitrarily placed rational point set."""
import json
from fractions import Fraction
from pathlib import Path
import sys

from flint import arb, ctx


def cross(a,b):
    return a[0]*b[1]-a[1]*b[0]


def rational(x):
    x=Fraction(x)
    return arb(x.numerator)/x.denominator


def verify(comp, pose, point_file, plan, output):
    sys.path.insert(0,str(comp/"inclusion-minimal-eleven"))
    from template import build
    meta=json.loads(plan.read_bytes())
    frozen=json.loads(pose.read_bytes())["center"]
    from math import ceil
    scale=Fraction(meta.get("translation_radius","0"))*2**44
    widths=[ceil(scale*(3 if i in (2,5,8) else 1)) for i in range(11)]
    box=[[x-z,x+z] for x,z in zip(frozen,widths)]
    poly,den=build(comp,False,meta["k"],True,box=box)
    poly=[tuple(map(int,q)) for q in poly]
    den=int(den)
    n=len(poly)
    assert meta["host_vertices"]==n and meta["den"]==den
    target=Fraction(meta["target"])
    twice=sum(cross(a,b) for a,b in zip(poly,poly[1:]+poly[:1]))
    allowance=2*den**2*target-twice
    assert allowance>0
    twice_poly=poly+poly
    prefix=[0]
    for a,b in zip(twice_poly,twice_poly[1:]):
        prefix.append(prefix[-1]+cross(a,b))
    candidates=set()
    for length in (1,2,3,4,6,8,12,16,22,30,40,52,68,86,110,142,180,230):
        if length>=n//2:
            continue
        for i in range(n):
            j=i+length;a,b=twice_poly[i],twice_poly[j]
            candidates.add((b[1]-a[1],a[0]-b[0],(allowance+prefix[j]-prefix[i])/den))
    rows=[(x,y,Fraction(b)) for x,y,b in meta["rows"]]
    assert all(w>0 and (x,y,w) in candidates for x,y,w in rows)
    points=json.loads(point_file.read_bytes())
    vertices=points["points"]; Q=points["denominator"]
    assert max((x-u)**2+(y-v)**2 for x,y in vertices for u,v in vertices)<=Q**2
    ctx.prec=192
    size=meta["grid"]
    intervals=[(2*arb.pi()*j/size).sin_cos() for j in range(size+1)]
    observed=set()
    minimum=arb(1)
    failures=[]
    for hand,j,ids,picks in meta["records"]:
        assert hand in (-1,1) and 0<=j<size and (hand,j) not in observed
        observed.add((hand,j))
        normals=[rows[i][:2] for i in ids]
        if len(normals)==2:
            a,b=normals
            assert cross(a,b)==0
            weights=[-a[0]*b[0]-a[1]*b[1],a[0]*a[0]+a[1]*a[1]]
        else:
            assert len(normals)==3
            a,b,c=normals
            weights=[cross(b,c),cross(c,a),cross(a,b)]
            if weights[0]<0: weights=[-x for x in weights]
        assert all(x>0 for x in weights)
        assert all(sum(w*vec[k] for w,vec in zip(weights,normals))==0 for k in [0,1])
        cs,ss=0,0; wall=Fraction(0); norm=arb(0)
        for i,k,w in zip(ids,picks,weights,strict=True):
            nx,ny,b=rows[i];x,y=vertices[k];y*=hand
            cs+=w*(nx*x+ny*y);ss+=w*(ny*x-nx*y);wall+=w*b
            norm+=w*arb(nx*nx+ny*ny).sqrt()
        ends=[((ss*s+cs*c)/Q-rational(wall))/norm for s,c in intervals[j:j+2]]
        gap=min(x.lower() for x in ends)
        if not gap>0:failures.append([hand,j,str(gap)])
        else: minimum=min(minimum,gap)
    assert len(observed)==2*size and size>=3
    result=dict(target=str(target),translation_radius=meta.get("translation_radius","0"),polygon_area=str(Fraction(twice,2*den**2)),host_vertices=n,point_count=len(vertices),physical_diameter_limit="1",cells=2*size,failed=len(failures),minimum_successful_clearance=str(minimum),first_failures=failures[:5])
    output.write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
    assert not failures
if __name__=="__main__":
    verify(*(Path(x).resolve() for x in sys.argv[1:]))
