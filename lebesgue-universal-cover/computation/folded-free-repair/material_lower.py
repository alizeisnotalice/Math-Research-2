"""Inner materials in the AH cap and a disjoint real face triangle for a specified host."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
from flint import arb,ctx,__version__ as flint_version


def run(comp, single, domain, mixed, difference):
    sys.path.insert(0,str(single));from certify_counterexample import geometry
    from certify_monotone_cost import posed_jets,quantize_chart,support_pins
    ns,body=geometry(comp)
    from certify_ah_tapered_minkowski_closure import geometry as get_geo
    from certify_nonmonotone_ah_core import hull,double_area,SCALE,LENGTH2_LO
    from certify_tapered_forced_support import rational as A
    old,_=get_geo(comp);ctx.prec=192;g=old['g']
    both=json.loads(mixed.read_bytes());instance=json.loads(domain.read_bytes())
    data=json.loads(difference.read_bytes());old_low=F(data['single']['lower'])
    cap_upper=F(both['cap_upper'])
    ta=g.plus(old['ns']['z'],old['ns']['D']);v=g.minus(old['q'](arb(1)),ta);xn=g.neg(g.rot(v))
    assert g.dot(v,v)>A(LENGTH2_LO)
    chart,outward,inward=quantize_chart(ns,g,ta,v,SCALE)
    tris=[(0,0),(SCALE,0)]; face=instance['extra_facet'];assert face in (3,4)
    normals=[]
    for i in range(-650,651):
        z=g.plus(xn,g.scale(arb(i)/1950,v));normals.append(g.scale(1/g.norm(z),z))
    def digest(p):return sha256(p.read_bytes()).hexdigest()
    files=[domain,mixed,difference,Path(__file__)]
    for name in ('pentagons.json','curves.json'):
        path=comp/'finite-b-core'/name; files.append(path)
        tab=json.loads(path.read_bytes())
        points=tab.get('combined_hull') or tab['vertices_ccw_integer'];den=tab.get('scale') or tab['grid']
        err=arb(0,arb(2)/den)
        for x,y in points:
            z=inward((arb(x)/den+err,arb(y)/den+err))
            if z is not None:tris.append(z)
    host=single/'single-H1-input.json';files.append(host);oldb=json.loads(host.read_bytes())['bodies']
    plane=ns['trig'](ns['A'](60*face));level=arb(0)
    for w,e,k,s,xx,zz in [z+['0','0'] for z in oldb]+instance['test_bodies']:
        vs,_=body(ns,F(w),F(e));es=[ns['sub'](vs[i],vs[j]) for i,j in ns['orders']]
        a=ns['A'](F(k));fn,t,q=support_pins(ns,vs,es,s,a)
        x,z=A(F(xx)),A(F(zz));t=g.plus(t,(x,(x+2*z)/arb(3).sqrt()))
        for nn in normals:
            j=inward(posed_jets(ns,vs,es,a,s,nn,t))
            if j is not None:tris.append(j)
        point=posed_jets(ns,vs,es,a,s,plane,t)
        if g.dot(xn,g.minus(point,ta))<0:
            level=max(level, (fn(plane)+g.dot(plane,t)-arb('.5')).lower())
    length=arb(both['old_supporting_segment_lengths'][face])
    assert length>arb('.41542') and level>0
    triangle_lb=F(int(float(level*length/2)*10**16)-1,10**16)
    assert A(triangle_lb)<length*level/2
    ah=F(double_area(hull(tris)),2*SCALE**2)*LENGTH2_LO
    lower=old_low+max(ah-cap_upper,F(0))+triangle_lb
    return dict(area_lower=str(lower),physical_AH=str(ah),old_AH_upper=str(cap_upper),extra_face=face,
                face_gain=str(triangle_lb),flint=flint_version,precision=ctx.prec,hashes={p.name:digest(p) for p in files})


if __name__=='__main__':
    p=ArgumentParser()
    for n in ('comp','single','instance','mixed','diff','out'):p.add_argument(n,type=Path)
    a=p.parse_args();result=run(a.comp.resolve(),a.single.resolve(),a.instance,a.mixed,a.diff)
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(float(F(result['area_lower'])))
