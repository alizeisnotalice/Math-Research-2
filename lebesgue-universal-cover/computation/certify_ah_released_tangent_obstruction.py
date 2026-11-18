"""Exclude every rigid pose for a larger triangular Minkowski perturbation of the pentagon."""
from fractions import Fraction as F
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from flint import arb,ctx
from certify_ah_released_tangent import load_geometry, triangle_support


def run(archive, long_patch):
    import audit_optimized_eh_chord as g
    ns=load_geometry(archive)
    A=ns['A']
    ctx.prec=160
    g.SIGMA=A('1.5497')*arb.pi()/180
    w=F(65,128) if long_patch else F('.42883')
    t=A('.0002') if long_patch else A('.000015')
    r=1/A(3).sqrt()
    h= A('.5')
    deg=arb.pi()/180
    v=ns['points'](A(w))[0]
    edges=[g.minus(v[k],v[l]) for k,l in ns['orders']]
    def pent(u):
        return ns['supports'](v,edges,(A(1),A(0)),u)
    def supp(u):
        return (1-t)*pent(u)+t*triangle_support(g,A,ns['reflection'](u))
    d=ns['direct']
    delta=A('.011') if long_patch else A('.018')
    pivot=g.plus(ns['z'],g.scale(A('.95' if long_patch else '.430'),ns['D']))
    new_normal=g.neg(g.rot(g.plus(g.scale(delta.cos(),d),g.scale(delta.sin(),g.rot(d)))))
    identity=supp(new_normal)-g.dot(new_normal,pivot)
    def tangent(w):
        x,dx=ns['point_jets'](A(w))
        pt=ns['reflection'](x[2])
        du=ns['reflection'](dx[2])
        n=g.scale(-1/g.norm(du),g.rot(du))
        return n,g.dot(n,pt)
    left,hl=tangent(F(54,128))
    right,hr=tangent(F(56,128))
    qa=supp(ns['reflection'](left))-hl
    if long_patch:
        offending,hshift=tangent(w)
        qa=supp(ns['reflection'](offending))-hshift-A('1e-6')
    assert identity>0 and qa>0
    corner=g.line_pair(left,right,hl,hr)
    hexagon,m,*rest=g.setup()
    caps=[m[3],m[5]]
    corners=[corner]
    if long_patch:
        rr=[tangent(F(53,128))]
        for k in range(57,125,4):
            n,offset=tangent(F(k,128))
            rr.append((n,offset+A('1e-6')))
        rr.append(tangent(F(1)))
        corners=[g.line_pair(n,m,a,b) for (n,a),(m,b) in zip(rr,rr[1:])]
    corner_P=min((h-g.dot(u,x)).lower() for x in corners for u in (*hexagon,*caps))
    assert corner_P>0
    other=[]
    for k in range(0,360,60):
        cos,sin=g.unit(-k*deg)
        for s in (1,-1):
            if (k,s) in [(0,1),(240,-1)]:
                continue
            def trans(u):
                return g.plus(g.scale(cos,(u[0],s*u[1])),g.scale(sin,g.rot((u[0],s*u[1]))))
            other.append(max((supp(trans(u))-h).lower() for u in caps))
    assert min(other)>0
    # Every integral number of sixty degrees is an exact compatibility root.
    # Uniform endpoint error |F(phi)-S sin(phi)| <= 3/2(1-cos(phi))
    # uses the three witnessed opposite contacts and the radius-one disks.
    slope=sum((-1)**k*g.dot(v[i],g.rot(hexagon[k])) for k,i in enumerate((1,3,0)))
    endmargin=-((1-t)*slope+3*r*t+A('1.5')*(A('.25')*deg).tan())
    assert endmargin>0
    worst=A(1)
    cells=0
    def evaluate(l,r,depth):
        nonlocal worst,cells
        mid=(l+r)/2
        z=A(mid)*deg
        span=A(l).union(A(r))*deg
        fm=sum((-1)**k*supp(g.unit(z+60*k*deg)) for k in range(3))-h
        derivative=sum((-1)**k*ns['supports'](v,edges,g.unit(span+60*k*deg),hexagon[0],True) for k in range(3))
        rad=A((r-l)/2)*deg
        err=(1-t)*rad*abs(derivative)+3*t*rad/A(3).sqrt()
        slack=(-fm-err).lower()
        if slack>0:
            worst=min(worst,slack)
            cells+=1
        else:
            assert depth<20,(l,r,slack)
            evaluate(l,mid,depth+1)
            evaluate(mid,r,depth+1)
    for i in range(236):
        l=F(1,2)+F(i,4)
        evaluate(l,l+F(1,4),0)
    return dict(w=str(w),long_patch=long_patch,rational_w=str(w),weight='0.0002' if long_patch else '0.000015',parameters_file='ah-released-tangent-long.json' if long_patch else 'ah-released-tangent-certificate.json',
                identity_violation=str(identity),QA_violation=str(qa),ten_other_pose_min_violations=str(min(other)),
                corner_P_clearance=str(corner_P),interior_angular_cells=cells,
                interior_negative_margin=str(worst),endpoint_slope_clearance=str(endmargin),
                complete_nonuniversality=True)


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('archive',type=Path)
    cli.add_argument('output',type=Path)
    cli.add_argument('--long',action='store_true')
    opt=cli.parse_args()
    path=opt.archive.resolve()
    sys.path.insert(0,str(path))
    result=run(path,opt.long)
    result['verifier_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    opt.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result)
