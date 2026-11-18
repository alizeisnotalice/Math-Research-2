"""Separate a moving curvature jump for arbitrarily many affinely independent true input bodies."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import arb, arb_series, ctx, __version__ as flint_version


def run(comp, certificate):
    sys.path.insert(0, str(comp))
    from deformed_reuleaux_support import body, environment
    base=json.loads(certificate.read_bytes())
    (wlo,whi),(elo,ehi)=[list(map(F,a)) for a in base['rectangle']]
    w=F('.81155')
    assert wlo<w<whi and elo<ehi
    ns=environment(str(comp),'1.5497');ctx.prec=192
    ebox=ns['box'](elo,ehi); body(ns,w,ebox)
    e=arb_series([ebox,1]);add,mul,sub,ref=map(ns.get,('add','mul','sub','reflection'))
    l,t=ns['chart'](add(ns['z'],mul(ns['A'](w),ns['D'])));p=ns['psi'](l);t=arb_series([t])-e
    vs=[add(ns['l0'],mul(l,ns['el'])),ref(add(ns['l0'],mul(p,ns['el']))),ref(ns['X'](p,-t)),
        add(ns['mid'],mul(t,ns['en'])),ns['X'](l,t)]
    vseries=[tuple(x if isinstance(x,arb_series) else arb_series([x]) for x in pair) for pair in vs]
    edge=[sub(vseries[i],vseries[j]) for i,j in ns['orders']]
    e0=[(x[0],y[0]) for x,y in edge]
    de=[(x[1],y[1]) for x,y in edge]
    dt=[ns['cross'](a,b) for a,b in zip(e0,de)]
    pick=next(j for j,d in enumerate(dt) if abs(d)>arb('.01'))
    speed=dt[pick];max_speed=abs(speed).upper()
    assert max_speed*ns['A'](ehi-elo)<arb.pi()
    sep= max(ns['dot'](e0[pick],v).upper() for j,v in enumerate(e0) if j!=pick)
    assert sep<1
    radii=[r for _,r in ns['sector']]
    assert len(radii)==len(edge) and all(r in (0,1) and r!=radii[j-1] for j,r in enumerate(radii))
    # The small residual budget also covers a uniform parallel neighborhood.
    # The rectangle body lies in the .58+1.4e-10 radius ball by the corner bound.
    assert 2*arb.pi()*(arb('.58')+arb('.00000000014'))<4
    pad=F('2.5e-9')
    cost=4*pad+4*pad**2
    assert cost<F(base['remaining_B'])
    return dict(w=str(w),e=[str(elo),str(ehi)],moving_knot=pick,edge_pair=ns['orders'][pick],
                angular_velocity=str(speed),other_knot_dot_upper=str(sep),radii=radii,
                neighborhood_radius=str(pad),neighborhood_area_upper=str(F(base['area_strict_upper'])+cost),
                rectangle_certificate_sha256=sha256(certificate.read_bytes()).hexdigest(),flint=flint_version,precision=ctx.prec,
                script_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('comp',type=Path);p.add_argument('rectangle',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    result=run(a.comp.resolve(),a.rectangle);a.output.write_text(json.dumps(result,indent=2)+'\n');print(result)
