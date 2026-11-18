"""A fixed triangular perturbation rules out deep AH cuts at all rotations and translations."""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import sys

import flint
from flint import ctx,arb
from certify_low_sigma_ah_budget import b,g,add,sub,mul,dot,cross,rot,neg,exact,pack,digest,model,pull,norm,CYCLE,OTHER,CASES
from certify_ah_released_tangent import triangle_support

W=F('.415')
DEPTH=F('.024')
T_WEIGHT=F('.002')
PIVOTS={'1.53':'.430295','1.52':'.431178','1.50':'.432974'}
PAIR_1={(0,1),(2,3),(4,5),(0,2),(0,6),(2,7),(4,6),(4,7)}

def verify(sigma,pivot):
    sigma=F(sigma)
    g.SIGMA=b.degree(sigma)
    half=g.HALF
    n,m,*_=g.setup()
    D=model.D_of(g.SIGMA)
    pts=b.anchors(D,exact(W))
    t=exact(T_WEIGHT)
    # The five exact sides and three opposite pairs have length one by construction.
    slack=min((1-dot(sub(pts[i],pts[j]),sub(pts[i],pts[j]))).lower()
              for i,j in combinations(range(len(pts)),2) if (i,j) not in PAIR_1)
    assert slack>0
    lens=b.lens(pts)
    R=1/arb(3).sqrt()
    vv=[(arb(0),-R),(-R*arb(3).sqrt()/2,R/2),(R*arb(3).sqrt()/2,R/2)]
    cones=[]
    for i in CYCLE:
        k=CYCLE.index(i)
        dl,dr=[sub(pts[i],pts[CYCLE[(k+s)%len(CYCLE)]]) for s in (-1,1)]
        if cross(dl,dr)<0:
            dl,dr=dr,dl
        assert cross(dl,dr)>0
        cones.append((pts[i],dl,dr))

    def base_lower(u):
        candidates=[dot(x,u).lower() for x in pts]
        for p,d,e in cones:
            if cross(d,neg(u))>=0 and cross(neg(u),e)>=0:
                candidates.append((1+dot(p,u)).lower())
        return max(candidates)

    def tri(u):
        return triangle_support(g,exact,neg(u))

    def lower(u):
        return (1-t)*base_lower(u)+t*tri(u).lower()

    def upper(u):
        return (1-t)*lens(u).upper()+t*tri(u).upper()

    endpoint=F('.1')
    slopes=[]
    for deg,sgn,active in [(0,1,(3,1,4)),(60,-1,(5,2,0))]:
        total=arb(0)
        for i in range(3):
            u=n[(deg//60+2*i)%6]
            point=pts[active[i]]
            assert (dot(u,point)-half).contains(0)
            # The triangle opposing face at each normal is an arc point or a
            # triangle tip. Use the exact support point, also at a fan seam.
            tip=max(vv,key=lambda z:dot(u,z).mid())
            far=min(vv,key=lambda z:dot(u,z).mid())
            if (dot(u,tip)-half).contains(0):
                vertex=tip
            else:
                assert (dot(u,far)+half).contains(0)
                vertex=add(far,u)
            assert (dot(u,vertex)-half).contains(0)
            total+=dot(mul(sgn,rot(u)),add(mul(1-t,point),mul(t,vertex)))
        slope=total* b.degree(endpoint).cos()-exact('1.5')*b.degree(endpoint).sin()
        assert slope>0,(sigma,'slope',deg,slope)
        slopes.append(slope.lower())
    margins=[]
    def scan(a,c):
        u=b.degree(a).union(b.degree(c))
        val=sum(lower(g.unit(u+2*i*g.PI/3)) for i in range(3))-exact('1.5')
        if val>0:
            margins.append(val.lower())
        else:
            assert c-a>F(1,1000000),(sigma,a,c,val)
            z=(a+c)/2
            scan(a,z)
            scan(z,c)
    scan(endpoint,F(60)-endpoint)
    # There are no internal compatibility roots. The perturbation still lies
    # within the two P in the dihedral orientations identity and Q_A.
    pwall=[]
    for k,r in ((0,1),(2,-1)):
        pwall.append(min((half-upper(pull(m[i],k,r))).lower() for i in (3,5)))
    assert min(pwall)>0,(sigma,pwall)
    dihedraleps=[]
    for k,r in OTHER:
        gap=max((lower(pull(m[i],k,r))-half).lower() for i in (3,5))
        assert gap>0,(sigma,k,r,gap)
        dihedraleps.append(gap)
    alpha=exact(DEPTH)
    ray=mul(-1/norm(D),D)
    cut=neg(rot(add(mul(alpha.cos(),ray),mul(alpha.sin(),rot(ray)))))
    X=add(b.z,mul(exact(pivot),D))
    ah=[]
    for k,r in ((0,1),(2,-1)):
        margin=lower(pull(cut,k,r))-dot(cut,X)
        assert margin>0,(sigma,k,r,margin)
        ah.append(margin.lower())
    return dict(sigma_degrees=str(sigma),ref_parameter=str(W),tilt=str(DEPTH),pivot=str(pivot),triangle_weight=str(T_WEIGHT),
                triangle_phase_degrees='60',diameter_slack_lower=pack(slack),end_slope_lower=[pack(r) for r in slopes],
                angular_cells=len(margins),minimum_angular_lower=pack(min(margins)),
                identity_and_QA_original_P_clearances=[pack(r) for r in pwall],
                other_ten_P_violations_lower=pack(min(dihedraleps)),
                identity_and_QA_AH_violations_lower=[pack(x) for x in ah],
                old_P_checks_use_lens_support_upper=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('archive',type=Path)
    ap.add_argument('output',type=Path)
    args=ap.parse_args()
    ctx.prec=192
    cases=[verify(s,F(PIVOTS[s])) for s in CASES]
    sources=['certify_moving_sigma_fold_splice.py','certify_two_sigma_ideal_placement.py','certify_ah_released_tangent.py']
    receipt=dict(complete=True,rows=cases,precision=ctx.prec,python=sys.version,flint=flint.__version__,
                 sources={f:digest(args.archive/f) for f in sources},
                 script_sha256=digest(Path(__file__)),
                 necessary_wall_source_sha256=digest(Path(__file__).with_name('certify_low_sigma_ah_budget.py')))
    args.output.write_text(json.dumps(receipt,indent=2)+'\n')
    for row in cases:
        print(row['sigma_degrees'],row['angular_cells'],
              [z['decimal'] for z in row['identity_and_QA_original_P_clearances']],
              [z['decimal'] for z in row['identity_and_QA_AH_violations_lower']],flush=True)
