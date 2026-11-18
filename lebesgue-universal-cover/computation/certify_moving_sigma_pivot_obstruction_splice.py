"""Uniform cap and tangent exclusions above the moving ideal compatibility fold."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
from flint import arb_series

import certify_moving_sigma_fold_splice as model
from certify_two_sigma_ideal_placement import add, anchors, cross, degree, dot, ea, g, mul, neg, qa, rot, sub

SIG,W,P=model.SIG,model.W,model.P
ANGLE_THRESHOLD=(F('1.54603709348333041715'),F('1.54603709348333041725'))
OTHER=tuple((k,s) for k in range(6) for s in (-1,1) if (k,s) not in ((0,1),(4,-1)))


def verify():
    sig_lo,sig_hi=SIG
    w_lo,w_hi=W
    deltas=tuple(model.delta(model.D_of(degree(s)),g.exact(P)) for s in ANGLE_THRESHOLD)
    assert deltas[0]>0 and deltas[1]<0
    nc=(64,256)
    cuts=[]
    projection=[]
    tangent=[]
    V=(g.arb(0),1/g.arb(3).sqrt())
    for i in range(nc[0]):
        s0,s1=[sig_lo+(sig_hi-sig_lo)*F(j,nc[0]) for j in (i,i+1)]
        s=degree(s0).union(degree(s1))
        sm=degree((s0+s1)/2)
        sd=ea((s0-s1)/2,(s1-s0)/2)
        mc=model.D_of(arb_series([sm,g.PI/180],2))
        dd=model.D_of(arb_series([s,g.PI/180,0],3))
        D=[m[0]+(m[1]+2*q[2]*sd)*sd for m,q in zip(mc,dd,strict=True)]
        shapes=[g.unit(g.PI*g.exact(F(1,6)+F(j,3))+s) for j in (3,5)]
        # The old m walls must exclude the same ten exact dihedral poses.
        directions={}
        for k,r in OTHER:
            c,ss=g.unit(-k*g.PI/3)
            def inverse(u):
                u=(u[0],r*u[1])
                return add(mul(c,u),mul(ss,rot(u)))
            directions[k,r]=list(map(inverse,shapes))
        for j in range(nc[1]):
            w0,w1=[w_lo+(w_hi-w_lo)*F(k,nc[1]) for k in (j,j+1)]
            w=ea(w0,w1)
            points=list(map(qa,anchors(D,w)))
            y=anchors(D,arb_series([w,1],2))[6]
            q,p=tuple(x[0] for x in y),tuple(x[1] for x in y)
            outward=mul(1/g.norm(p),neg(rot(p)))
            projection.append(dot(outward,sub(V,q)).lower())
            tangent.append(g.norm(p).lower())
            cuts.append(min(max((dot(u,x)-g.HALF).lower() for u in d for x in points)
                            for d in directions.values()))
    assert min(min(cuts),min(tangent),min(projection))>0
    return {'complete':True,'sigma_degrees':list(map(str,SIG)),'w':list(map(str,W)),
            'pivot_parameter':str(P),'pivot_crossing_sigma_degrees':list(map(str,ANGLE_THRESHOLD)),
            'pivot_crossing_Delta_signs':[str(d) for d in deltas],
            'cap_and_tangent_parameter_cells':list(nc),
            'strict_other_ten_cap_violation_lower':str(min(cuts)),
            'strict_triangle_projected_excess_lower':str(min(projection)),
            'tangent_norm_lower':str(min(tangent)),
            'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}

if __name__=='__main__':
    Path(sys.argv[1]).write_text(json.dumps(verify(),indent=2)+'\n')
