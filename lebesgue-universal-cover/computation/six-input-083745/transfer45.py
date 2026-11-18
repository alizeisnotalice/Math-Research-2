"""Audit all geometry for transferring the 0.83745 original-tree target to fixed rational sets."""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import argparse
import json
import sys
import tarfile


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ('comp', 'tree', 'table32', 'table64', 'scope_tar', 'out'):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    assert __debug__
    sys.path.insert(0,str(a.comp.resolve()))
    from certify_polygon_083731 import ORDERS, ar, arb, ctx, lower_radius
    from certify_tomar_six_high import SourceTree
    from certify_refined_material_contraction import verify
    from certify_shallow_cap_transfer import diameter_squared
    from tangent_polygons_finite_cover import generate
    from select_tomar_083731_scores import radius_rows
    from tomar_d3_mask import CUT
    ctx.prec=160
    read=lambda name: json.loads((a.comp/name).read_bytes())
    goal=F(120689264254125501,2**57)
    scale,r,R=F('.999999006'),F(1,2),F('.703')
    rho,q,h,wall,D=map(F, ('.000000860','.528619','.0125','.47184','.495'))
    verified=verify(a.comp,a.tree,a.table32,a.table64)
    assert verified==read('refined-material-contraction.json') and F(verified['curve_scale'])==scale
    tree=SourceTree(a.tree)
    with tarfile.open(a.scope_tar,'r:gz') as t:
        subset=json.load(t.extractfile('inputs/ordinary/scope.json'))
    keys=sorted(tuple(z) for z in subset['leaves']+subset['deferred'])
    assert len(set(keys))==len(keys)==3273
    maxima=[(sum(max(abs(lo),abs(hi))**2 for lo,hi in tree.leaf(*k)[:2]),k) for k in keys]
    max_square,key=max(maxima)
    assert max_square<(r-wall)**2
    pi=arb.pi(); sr=scale*r
    assert all(lower_radius(n,pi)>ar(wall) for n in (5,7,9))
    disk_point=(pi-ar(sr/R).acos())*ar(sr)**2+ar(sr)*ar(R**2-sr**2).sqrt()
    assert disk_point>ar(goal) and R**2<2*sr**2
    root=tree.root
    assert root['ns']==[3,5,7,9] and root['split0']==[2]*11
    assert all(F(root['lo'][j])==0 for j in (2,5,8))
    assert pi/15<ar(CUT)<ar(F(root['hi'][2]))
    assert all(ar(F(root['hi'][j]))>2*pi/n for j,n in ((5,7),(8,9)))
    exterior=[]
    seed=read('tomar-083745-future-seed.json'); radii=radius_rows(read('five-body-radial-affine.json'),F('.83745'))
    assert seed['radii']==radii and seed['target_integer']==goal*2**57
    for i,n in enumerate(ORDERS[1:]):
        j=3*i
        displacement=min(F(root['hi'][j]),F(root['hi'][j+1]),-F(root['lo'][j]),-F(root['lo'][j+1]))
        inside=ar(scale)*lower_radius(n,pi)
        area=pi*(ar(sr)**2+inside**2)/2+ar(displacement)*(ar(sr)+inside)
        assert area>ar(goal) and F(radii[i]['integer_area_lower'],2**57)>goal
        assert F(radii[i]['integer_radius'],2**44)<min(displacement,F('.069'))
        exterior.append(str(area))
    delta=F('.6')*(1-scale);loss=F(22,7)*(2*R*delta+delta**2)
    assert F('.837455')-loss>goal
    five=read('five-global-minimum-single-worker-replay.json')
    assert five['complete'] and arb(five['minimum_interval'])>ar(F('.837295'))
    assert D*D>R*R-sr*sr
    cap=scale**2*F('.837295')+sr*h*h/D-goal
    assert cap>0
    cloud=read('finite-point-leaf-input.json'); den=cloud['denominator']; assert len(cloud['points'])==30
    assert all((F(x,den)+F('.01007'))**2+F(y,den)**2<q*q for x,y in cloud['points'])
    m=rho*q/(1-scale);v=(1+scale)*r*r+(1-scale)*m*m
    slope=1+scale-(2*m*r+v)/(4*r*r)
    margin=4*(wall-m)*scale*scale-(1+scale)*h-max(0,(v-2*m*(r+h))/(2*r+h))
    assert m<wall<r and slope>0 and margin>0
    dsq=diameter_squared(cloud); t=F(isqrt(10**36*dsq.denominator//dsq.numerator),10**18)
    factor=t/(1-rho)
    groups,budget=generate((scale*factor)**2,F(1))
    assert F(budget['common_dilation'])*scale*factor<1
    groups[-1]={'body':'P30','denominator':t.denominator*den,'points':[[t.numerator*x,t.numerator*y] for x,y in cloud['points']]}
    assert diameter_squared(groups[-1])<=1
    raw=(json.dumps(groups,separators=(',',':'))+'\n').encode()
    assert raw==(a.comp/'native-wall-083744-points.json').read_bytes()
    assert sum(len(g['points']) for g in groups)==17513
    out=dict(status='geometry for transfer at G; requires validated area trees and full-domain seed',
        goal=str(goal),radial_numerators=[x['integer_area_lower'] for x in radii],root_exterior=exterior,disk_point_outer=str(disk_point),
        scale=str(scale),point_shrink=str(rho),height=str(h),wall=str(wall),attaining_leaf=key,native_leaves=len(keys),translation_squared=str(max_square),
        cap_height_margin=str(cap),cap_slope=str(slope),cap_comparison_margin=str(margin),pure_steiner_margin=str(F('.837455')-loss-goal),material_squared_maxima=verified['material_squared_maxima'],
        coordinates_sha256=sha256(raw).hexdigest(),factor=str(factor),endpoint=str(factor**2*goal),display=float(factor**2*goal),point_budget=budget)
    a.out.write_text(json.dumps(out,indent=2)+'\n')
    print(out['endpoint'], out['display'], flush=True)

if __name__=='__main__':
    main()
