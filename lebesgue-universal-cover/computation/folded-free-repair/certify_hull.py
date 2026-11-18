"""Certify a fixed-placement enlargement of C* by one exact common convex hull."""
import argparse
import json
from pathlib import Path
import platform
from fractions import Fraction as F
import sys
from flint import arb,ctx,__version__ as flint_version
from certify_counterexample import digest,geometry
from certify_joint import dependency_log
from certify_monotone_cost import posed_jets,support_pins

def prove(comp,host,pair,single_repair,comparison):
    sys.path.insert(0,str(comp))
    sys.path.insert(0,str(single_repair))
    sys.path.insert(0,str(comparison))
    manifest=dependency_log({'computation':comp,'single-repair':single_repair,'compare':comparison})
    import certify_single_repair as repair
    from certify_ah_tapered_mirror_counterexample import Support,Polygon
    from certify_tapered_forced_support import Profile
    from compare import materialize
    from certify_nonmonotone_ah_core import hull,double_area
    from certify_tapered_forced_support import rational as A
    ns,body=geometry(comp)
    old,core,container=repair.refined_core(comp,4096,hull)
    ctx.prec=192
    p=json.loads(pair.read_bytes())
    h=json.loads(host.read_bytes())
    folder=comp/'finite-b-core'
    pent=json.loads((folder/'pentagons.json').read_bytes())
    curves=json.loads((folder/'curves.json').read_bytes())
    n=pent['directions']
    D=pent['grid']
    Q=repair.Q
    assert n==curves['N'] and D==curves['scale'] and Q%D==0
    assert F(pent['coordinate_rounding_error_strict'])<=F(2,D)
    assert F(curves['coordinate_error_asserted'])<=F(2,D)
    sx=Q//D
    coarse=hull([(sx*x,sx*y)for x,y in pent['vertices_ccw_integer']+curves['combined_hull']])
    assert h['translation_rule']=='pin_hexagon_faces_0_and_2'
    first,last=[h['family_interval'][key] for key in ('first','last')]
    bones=h['bodies']
    assert len(bones)==3 and (first,last)==(1,2)
    assert bones[first][0]==bones[last][0] and bones[first][2:]==bones[last][2:]
    R=A(F(12,8))*(A(F(bones[last][1])-F(bones[first][1])))**2
    assert R<arb('8.86e-12')
    objects=[]
    displacements=[]
    release=[F(0),F(p.get('eta','.000000301')),F(0),F(0),F(0),F(0)]
    extra=p.get('extra_facet')
    if extra is not None:
        assert extra in (3,4)
        release[extra]=F('.000000020') if extra==3 else F('.000000321')
    for j,recipe in enumerate(bones+p['test_bodies']):
        w,e,d,s,*shift=recipe
        assert s in (1,-1) and len(shift) in (0,2)
        vs,_=body(ns,F(w),F(e))
        es=[ns['sub'](vs[i],vs[k])for i,k in ns['orders']]
        turn=ns['A'](F(d))
        h0,t,q=support_pins(ns,vs,es,s,turn)
        if j>=3:
            x,z=(A(F(x)) for x in shift) if shift else (arb(0),arb(0))
            t=(t[0]+x,t[1]+(x+2*z)/arb(3).sqrt())
            heights=(x,q+x+z,z,-x,-q-x-z,-z)
            assert all(gap<=A(bound) for gap,bound in zip(heights,release)),(j,list(map(str,heights)))
            if extra is None:
                assert 0<q<A(F(p['eta'])),(j,str(q))
            displacements.append(list(map(str,heights)))
        infl=R if j in (first,last) else arb(0)
        objects.append((vs,es,s,turn,h0,t,q,infl))
    v0=old['g']
    new=[]
    m=1440
    for j in range(m):
        normal=ns['trig'](ns['A'](F(360*j,m)))
        for vs,es,s,turn,_,t,_,_ in objects:
            x,y=posed_jets(ns,vs,es,turn,s,normal,t)
            new.append((round(float(x.mid())*Q),round(float(y.mid())*Q)))
    seed=hull(core+coarse+tuple(new))
    norms=[(v-y,x-u)for (x,y),(u,v) in zip(seed,seed[1:]+seed[:1])]
    old_core,old_pent,container_bound=[repair.support_scan(V,norms) for V in (core,coarse,container)]
    placements=json.loads((folder/'host-poses.json').read_bytes())['poses']
    datasets={k:json.loads((comp/k).read_bytes())for k in {r['source']for r in placements}}
    mirror=Support(datasets['ah-tapered-mirror-counterexample-input.json'],arb.pi())
    mirrors={'Support':mirror,'Polygon':Polygon(mirror)}
    helpers=[]
    for row in placements:
        if 'w' in row:
            _,fn=body(ns,F(row['w']),F(row['e']))
        elif 'body_type' in row:
            fn=mirrors[row['body_type']].support_direction
        else:
            fn=Profile(datasets[row['source']][row['body_index']]).support
        cos,sin=ns['trig'](ns['A'](F(row['turn'])))
        helpers.append((cos,sin,row['reflection'],fn))
    pad=2*(arb.pi()/(2*n)).sin()**2+arb(3)/D
    oblique=[old['tilt'][j] for j in (3,5)]
    oblique_cost=[]
    for vs,es,s,turn,hh,t,q,_ in objects[3:]:
        oblique_cost.append([str(hh(U)+v0.dot(U,t)-arb('.5')) for U in oblique])
    assert all(arb(gap)<0 for row in oblique_cost for gap in row)
    bounds=[]; refinements=0
    for (u,v),core_h,polygon_h,container_h in zip(norms,old_core,old_pent,container_bound):
        d=arb(u*u+v*v).sqrt()
        U=(arb(u)/d,arb(v)/d)
        old_h=min(arb(container_h),max(arb(core_h),(arb(polygon_h)+pad*d*Q).upper()))
        if old_h>core_h:
            ux,uy=U
            try:
                exact=max(hh((co*ux+s*si*uy,-si*ux+s*co*uy)).upper()
                          for co,si,s,hh in helpers)
                old_h=min(old_h,max(arb(core_h),(exact*d*Q).upper()))
                refinements+=1
            except AssertionError:
                pass
        h_new=max((hh(U)+v0.dot(U,t)+I).upper() for _,_,_,_,hh,t,_,I in objects)
        raw=max(old_h,(h_new*d*Q).upper())
        bound=int(float(raw))+5
        while not arb(bound)>raw:
            bound+=8
        bounds.append(bound)
    outer=materialize(norms,bounds,hull,repair.outer_rational)
    up=F(double_area(outer),2*Q**2)
    assert up<F(p['absolute_target']),(float(up),float(F(p['absolute_target'])))
    scale=10**22
    ceiling=up.numerator*scale//up.denominator+1
    assert up<F(ceiling,scale)
    decimal=f'{ceiling//scale}.{ceiling%scale:022d}'
    return dict(complete=True,area_upper=str(up),area_upper_decimal=decimal,exact_old_refinements=refinements,
                all_repository_inputs=manifest(),six_displacements=displacements,facet_budget=list(map(str,release)),
                slack=str(F(p['absolute_target'])-up),angles=[b[2:]for b in p['test_bodies']],
                width_excess=[str(z[6])for z in objects[3:]],oblique_excess=oblique_cost,
                halfplanes=len(norms),outer_vertices=len(outer),old_sample_count=n,old_sample_radius=str(pad),
                input_sha256={str(x.name):digest(x)for x in (
                   host,pair,single_repair/'certify_single_repair.py',comparison/'compare.py',folder/'pentagons.json',folder/'curves.json',
                   Path(__file__),comparison/'certify_monotone_cost.py',comparison/'certify_counterexample.py')},
                python=platform.python_version(),flint=flint_version,precision=ctx.prec)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for f in ('comp','host','pair','single_repair','comparison','out'):ap.add_argument(f,type=Path)
    a=ap.parse_args()
    a.out.write_text(json.dumps(prove(a.comp.resolve(),a.host,a.pair,a.single_repair.resolve(),a.comparison.resolve()),indent=2)+'\n')
