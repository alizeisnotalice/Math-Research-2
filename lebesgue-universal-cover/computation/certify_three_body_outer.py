"""Exact rational Steiner outer certificate for one complete D,B3,B5 placement."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import platform

import flint
from flint import arb,ctx
from certify_rotation_inballs import hull

GRID=10**12


def certify(params,samples):
    assert len(params)==5 and 64<=samples<=131072
    ctx.prec=112;pi=arb.pi();params=list(map(Q,params))
    def ball(q):return arb(q.numerator)/q.denominator
    def integer(x):
        n=round(float(x.mid())*GRID)
        assert abs(x*GRID-n)<1
        return n
    p=list(map(ball,params)); points=[]
    def emit(center,radius,first,step,count):
        for i in range(count):
            sn,cs=(first+i*step).sin_cos()
            points.append((integer(center[0]+radius*cs),integer(center[1]+radius*sn)))
    emit((0,0),arb(1)/2,arb(0),pi/(4*samples),8*samples)
    for n,rho,x,y in [(3,arb(0),p[0],p[1]),(5,p[2],p[3],p[4])]:
        R=1/(2*(pi/(2*n)).cos())
        for j in range(n):
            angle=rho+2*pi*j/n;sn,cs=angle.sin_cos()
            center=(x+R*cs,y+R*sn)
            emit(center,arb(1),angle+pi-pi/(2*n),pi/(n*samples),samples+1)
    polygon=hull(points);assert len(polygon)>=3
    twice=0;per=0
    for (a,b),(u,v) in zip(polygon,polygon[1:]+polygon[:1]):
        twice+=a*v-b*u
        square=(u-a)**2+(v-b)**2;ceil=isqrt(square);ceil+=ceil**2<square
        assert ceil**2>=square;per+=ceil
    area=Q(twice,2*GRID**2);perimeter=Q(per,GRID)
    # Arc interpolation error <= r*(angle step)^2/8; endpoint vector error <2/GRID.
    pi_hi=Q(22,7);assert pi<ball(pi_hi)
    eta=Q(2,GRID);sagitta=pi_hi**2/(72*samples**2);epsilon=eta+sagitta
    upper=area+perimeter*epsilon+pi_hi*epsilon**2
    # Conversely polygon lies in H+eta*B <= (1+2*eta)*H because H contains D.
    lower=area/(1+2*eta)**2
    assert lower<upper
    digest=sha256()
    for u,v in polygon:digest.update(f'{u},{v}\n'.encode())
    return {'python':platform.python_version(),'python_flint':flint.__version__,'precision_bits':ctx.prec,
            'parameters':[str(x) for x in params],'K':samples,'grid':GRID,'samples':len(points),'hull_vertices':len(polygon),
            'ordered_hull_sha256':digest.hexdigest(),'polygon_area':str(area),'perimeter_upper':str(perimeter),'sagitta_upper':str(sagitta),
            'total_outer_offset':str(epsilon),'placement_area_lower':str(lower),'placement_area_upper':str(upper),
            'interval_display':[float(lower),float(upper)],'width_display':float(upper-lower),
            'gap_above_0.83475':str(upper-Q(3339,4000)),
            'scope':'two-sided area bound for this exact decimal placement; its upper endpoint also bounds M(D,B3,B5)'}

if __name__=='__main__':
    ap=ArgumentParser();ap.add_argument('probe',type=Path);ap.add_argument('--K',type=int,default=8192);a=ap.parse_args()
    params=json.loads(a.probe.read_bytes())['recommended_decimal_parameters'];print(json.dumps(certify(params,a.K),indent=2))
