"""Build diameter-one rational groups for the .83742001 target from the .83742 contracted proof."""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("computation","shallow","coordinates","receipt"):
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    assert __debug__ and not args.coordinates.exists() and not args.receipt.exists()
    sys.path.insert(0,str(args.computation.resolve()))
    from certify_finite_083740_transfer import certified_core_scale
    from certify_polygon_083731 import ar, arb, ctx, diameter_squared
    from tangent_polygons_finite_cover import tangent_points
    ctx.prec=160
    shallow=json.loads(args.shallow.read_bytes())
    assert shallow['complete'] and Fraction(shallow['decimal_floor'])>Fraction('.83742')
    assert shallow['high_margin_sha256']=='b3632c4234469d95626091c56e96adfa9ee501c9237f6537bc86ee6e6473aa62'
    G=Fraction(shallow['integer_goal'])
    assert G==Fraction(120684940798483225,2**57)
    s=Fraction('.999999035');rho=Fraction('.000000735');T=Fraction('.83742001')
    f=Fraction(shallow['rational_scale'])
    groups=[];fits=[];diameters=[]
    orders=(0,3,5,7,9);samples=(1965,379,227,163,127)
    for n,m in zip(orders,samples,strict=True):
        points=tangent_points(n,m)
        assert len(points)==(n*(m+1) if n else m)
        maximum=max((x-a)**2+(y-b)**2 for (x,y),(a,b) in combinations(points,2))
        den=isqrt(maximum-1)+1
        assert (den-1)**2<maximum<=den**2
        kappa,interval,edges=certified_core_scale(n,points,den)
        f=min(f,kappa/s)
        fits.append(dict(body=n,samples=m,edges=edges,largest_squared_distance=maximum,
                         denominator=den,core_scale=str(kappa),arb_support_lower=interval))
        groups.append(dict(body=n,denominator=den,points=points))
        diameters.append(Fraction(maximum,den**2))
    assert f*f*G>T and all(f*s<=Fraction(z['core_scale']) for z in fits)
    cloud=json.loads((args.computation/'finite-point-leaf-input.json').read_bytes())
    scale=f*(1-rho)
    groups.append(dict(body='P30',denominator=cloud['denominator']*scale.denominator,
                       points=[[scale.numerator*x,scale.numerator*y] for x,y in cloud['points']]))
    diameters.append(diameter_squared(groups[-1]))
    assert max(diameters)<=1
    raw=(json.dumps(groups,separators=(',',':'))+'\n').encode()
    args.coordinates.write_bytes(raw)

    pi=arb.pi()
    dilation=1/(ar(s)*ar(T/G).sqrt())
    disk=next(k for k in range(3,2500) if (pi/(2*k)).cos()/(pi/k).cos()<dilation)
    assert (pi/(2*(disk-1))).cos()/(pi/(disk-1)).cos()>dilation
    lower_counts=[disk]
    for n in orders[1:]:
        exact=pi/(2*n*(1/dilation).acos())
        k=int(exact.lower().floor().unique_fmpz())
        assert k<exact<k+1
        lower_counts.append(n*k)
    lower_counts.append(30)
    sizes=[len(g['points']) for g in groups]
    result=dict(complete=True,points_sha256=sha256(raw).hexdigest(),groups=sizes,total=sum(sizes),
                contracted_goal=str(G),common_factor=str(f),sixth_scale=str(scale),
                exact_area_floor=str(G*f*f),decimal_floor='.83742001',fits=fits,
                squared_diameters=list(map(str,diameters)),containment_method_vertex_lower=lower_counts,
                method_total_lower=sum(lower_counts),shallow_receipt_sha256=sha256(args.shallow.read_bytes()).hexdigest())
    args.receipt.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('total','groups','exact_area_floor','method_total_lower','points_sha256')},indent=2))


if __name__=='__main__':
    main()
