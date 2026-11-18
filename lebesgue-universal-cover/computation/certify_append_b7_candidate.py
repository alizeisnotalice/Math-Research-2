"""Feasible B7 extension and Hausdorff stability at the exact three-body minimizer."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import platform

import flint
from flint import arb,ctx
from certify_append_b7_at_three_minimum import A,B
from certify_rotation_inballs import hull
from certify_three_body_symmetric_slice import interval,q,triple

GRID=10**12
B7_PLACEMENT=tuple(map(F,('0.2512845599997','-0.0003342551022','0.0126730937674')))
SNAP=F(1,10**12)
TAU_NOMINAL=F('0.837130000010')
TARGET=F('0.8363')


def close_exact_reference():
    a=interval('-.012154558715405409','-.012154558715405408')
    b=triple(a)[1]
    assert abs(a-q(A))<q(SNAP) and abs(b-q(B))<q(SNAP)
    assert TAU_NOMINAL/(1+2*SNAP)**2 > F('0.83713')
    assert 2*F(99,70)>0 and F(99,70)**2>2
    cos18=(arb.pi()/10).cos()
    assert 1/(2*cos18)<q('.526')
    lipschitz=F('1.941')
    assert F(99,70)+F('.526')<lipschitz
    candidate=F(1,7870)
    eta=SNAP+lipschitz*candidate
    low=TAU_NOMINAL/(1+2*eta)**2
    assert low>TARGET
    return a,b,candidate,low,lipschitz


def feasible_area(samples):
    pi=arb.pi(); pts=[]
    def integer(v):
        n=round(float(v.mid())*GRID)
        assert abs(v*GRID-n)<1
        return n
    def emit(center,r,angle,step,count):
        for j in range(count):
            s,c=(angle+j*step).sin_cos()
            pts.append((integer(center[0]+r*c),integer(center[1]+r*s)))
    emit((0,0),arb(1)/2,arb(0),pi/(4*samples),8*samples)
    w,x,y=map(q,B7_PLACEMENT)
    for n,a,center in [(3,arb(0),(q(A),0)),(5,arb(0),(q(B),0)),(7,w,(x,y))]:
        radius=1/(2*(pi/(2*n)).cos())
        for j in range(n):
            ang=a+2*pi*j/n
            s,c=ang.sin_cos();v=(center[0]+radius*c,center[1]+radius*s)
            emit(v,arb(1),ang+pi-pi/(2*n),pi/(n*samples),samples+1)
    pg=hull(pts);twice=0;perim=0
    for p,v in zip(pg,pg[1:]+pg[:1]):
        twice+=p[0]*v[1]-p[1]*v[0]
        d=(p[0]-v[0])**2+(p[1]-v[1])**2;l=isqrt(d);l+=l*l<d
        assert l*l>=d
        perim+=l
    area=F(twice,2*GRID**2);perimeter=F(perim,GRID)
    pi_hi=F(22,7);assert pi<q(pi_hi)
    outer_offset=F(2,GRID)+pi_hi**2/(72*samples**2)+SNAP
    outer=area+perimeter*outer_offset+pi_hi*outer_offset**2
    assert outer<F('0.837133895')
    digest=sha256()
    for x,y in pg: digest.update(f'{x},{y}\n'.encode())
    return {'sampling_steps':samples,'samples':len(pts),'vertices':len(pg),'grid':GRID,
            'hull_sha256':digest.hexdigest(),'area':str(area),'perimeter_upper':str(perimeter),
            'offset':str(outer_offset),'true_reference_placement_area_upper':str(outer)}


def main():
    parser=ArgumentParser();parser.add_argument('--samples',type=int,default=16384);args=parser.parse_args()
    assert args.samples>=64
    a,b,width,low,c=close_exact_reference()
    exact_lo,exact_hi=map(F,('.8347809457752035','.8347809457752037'))
    area0lo,area0hi=F('.83713'),F('.837133895')
    gainlo,gainhi=area0lo-exact_hi,area0hi-exact_lo
    print(json.dumps({'python':platform.python_version(),'flint':flint.__version__,'precision_bits':ctx.prec,
       'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'nominal_base':[str(A),'0','0',str(B),'0'],
       'exact_three_min_a_b':[str(a),str(b)],'base_snap_error_bound':str(SNAP),
       'required_nominal_lower':str(TAU_NOMINAL),'test_B7_angle_tx_ty':list(map(str,B7_PLACEMENT)),
       'exact_three_minimum_plus_B7_bracket':[str(area0lo),str(area0hi)],'necessary_increment_bracket':[str(gainlo),str(gainhi)],
       'coordinate_halfwidth_around_exact_minimum':str(width),'hausdorff_constant':str(c),
       'four_body_area_on_box_strictly_above':str(TARGET),'certified_moving_box_lower':str(low),
       'outer_area':feasible_area(args.samples)},indent=2))


if __name__=='__main__':
    main()
