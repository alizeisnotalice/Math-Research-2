"""Put a large original-coordinate box inside a certified contact-coordinate basin."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform

import flint
from flint import arb,ctx
from four_contact_taylor import Model
from certify_four_body_basin import CENTERS,RAW_CENTER
from certify_four_body_contact_stationary import add,geometry,pi,q,r,regular,unit_root


def certify(widths,radius,transverse,degree):
    Model.degree=degree
    x=Model({(0,0):q(RAW_CENTER[0]),(1,0):q(widths[0])})
    y=Model({(0,0):q(RAW_CENTER[1]),(0,1):q(widths[1])})
    na=unit_root(add((x,y),regular(3,0)),r,1)
    nb=unit_root(add((x,y),regular(3,1)),r,-1)
    a=(na[1]/na[0]).atan();b=pi/2-(nb[0]/nb[1]).atan()
    real=arb(0,1)
    assert na[0].evaluate(real,real)>0 and nb[1].evaluate(real,real)>0
    offsets=[(v-q(c)).evaluate(real,real).abs_upper() for v,c in zip((a,b),CENTERS)]
    assert max(offsets)<q(radius),offsets
    data=geometry(a,b)[0][2:]
    shifts=[(m-q(c)).evaluate(real,real).abs_upper()+q(w) for m,c,w in zip(data,RAW_CENTER[2:],widths[2:])]
    assert max(shifts)<q(transverse),shifts
    return offsets,shifts


def main():
    parser=ArgumentParser()
    parser.add_argument('basin',type=Path)
    parser.add_argument('--width',type=F,nargs='+',default=[F('0.000005')])
    parser.add_argument('--degree',type=int,default=10)
    args=parser.parse_args()
    widths=args.width*8 if len(args.width)==1 else args.width
    assert len(widths)==8 and min(widths)>0
    raw=args.basin.read_bytes(); basin=json.loads(raw)
    assert list(CENTERS)==basin['tangent_center'] and list(map(str,RAW_CENTER))==basin['raw_center']
    source=Path(__file__)
    assert basin['script_sha256']==sha256(source.with_name('certify_four_body_basin.py').read_bytes()).hexdigest()
    radius=F(basin['tangent_halfwidth']); transverse=F(basin['transverse_halfwidth'])
    a,b=certify(widths,radius,transverse,args.degree)
    print(json.dumps({'python':platform.python_version(),'flint':flint.__version__,'precision_bits':ctx.prec,
                      'script_sha256':sha256(source.read_bytes()).hexdigest(),'basin_sha256':sha256(raw).hexdigest(),
                      'taylor_degree':args.degree,'raw_center':list(map(str,RAW_CENTER)),'raw_halfwidths':list(map(str,widths)),
                      'tangent_offset_upper':list(map(str,a)),'transverse_offset_upper':list(map(str,b))},indent=2))


if __name__=='__main__':
    main()
