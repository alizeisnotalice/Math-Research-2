"""Reject the entire H, two tilted caps and one deep AH half-plane intersection."""
from argparse import ArgumentParser
from hashlib import sha256
import json
from pathlib import Path
import sys
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_tapered_forced_support import Profile,rational


def run(comp,profile):
    assert sha256((HERE/'certify_tapered_forced_support.py').read_bytes()).hexdigest()=='2cb253f3eddf195dbb05c65a84d9a9c40d17e8101f8962348f1e27e80de6d646'
    sys.path.insert(0,str(comp.resolve()))
    from certify_ah_tapered_minkowski_closure import geometry
    old,_=geometry(comp);g,ns=old['g'],old['ns']
    ctx.prec=192
    p,delta=rational('.8697'),rational('.01262')
    angle=delta
    normal=g.neg(g.rot(g.plus(g.scale(angle.cos(),ns['direct']),g.scale(angle.sin(),g.rot(ns['direct'])))))
    W=g.plus(ns['z'],g.scale(p,ns['D']))
    walls={k:(u,h)for k,u,h in old['active']if k in ('K-C3','I-E3')}
    walls['deep']=(normal,g.dot(normal,W))
    labels=['K-C3','deep','I-E3','I-E3','I-E3','I-E3','I-E3','I-E3','deep','K-C3','K-C3','K-C3']
    spec,=json.loads(profile.read_bytes())
    f=Profile(spec);mid,ends=f.root_clearance()
    rows=[];deep=[];caps=[]
    for pose,label in zip(((k,s)for k in range(0,360,60)for s in (-1,1)),labels,strict=True):
        u,h=walls[label];v=f.support(old['preimage'](u,*pose))-h
        assert v>rational('0.00000857')
        rows.append(dict(pose=pose,wall=label,violation=str(v)))
        (deep if label=='deep' else caps).append(v.lower())
    # A fixed direction, rotated AH edge translates linearly with p.
    coeff=g.dot(normal,ns['D'])
    assert not coeff.contains(0)
    pthreshold=p+min(deep)/coeff
    return dict(complete=True,p='.8697',delta='.01262',poses=rows,min_deep=str(min(deep)),min_cap=str(min(caps)),
                pivot_coefficient=str(coeff),pivot_at_first_escape=str(pthreshold),root_interior_clearance=mid,
                endpoint_slope_clearance=ends,profile_sha256=sha256(profile.read_bytes()).hexdigest(),
                parameters=list(map(str,f.x)),det=str(f.det))


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    for n in ('computation','profile','output'):parser.add_argument(n,type=Path)
    args=parser.parse_args()
    if not __debug__:raise SystemExit('Assertions must be enabled.')
    out=run(args.computation,args.profile)
    out['verifier_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(out['min_deep'],out['min_cap'],out['pivot_coefficient'],out['pivot_at_first_escape'])
