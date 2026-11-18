"""Exclude a two-parameter family of genuine Reuleaux pentagons from the posed B core."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
import flint
from flint import arb,ctx


def run(comp):
    comp=comp.resolve()
    sys.path.insert(0,str(comp))
    from certify_ah_tapered_minkowski_closure import geometry
    from deformed_reuleaux_support import body,environment
    from certify_tapered_forced_support import Profile
    from certify_ah_tapered_mirror_counterexample import Support,Polygon
    from certify_nonideal_pentagon_pose_library import verify
    old,_=geometry(comp)
    ns=environment(str(comp),'1.5497')
    ctx.prec=160
    g=old['g']
    A=ns['A']
    ta=g.plus(ns['z'],ns['D'])
    tp=old['q'](arb(1))
    chord=g.minus(tp,ta)
    chord_normal=g.neg(g.rot(chord))
    line_normal=g.neg(g.rot(ns['direct']))
    assert g.cross(line_normal,chord_normal)>0
    assert g.dot(chord_normal,g.neg(ta))<0
    angles=['112.325','113.311']
    directions=[ns['trig'](A(F(s))) for s in angles]
    for n in directions:
        assert g.cross(line_normal,n)>0 and g.cross(n,chord_normal)>0
    caps=[old['tilt'][i] for i in (3,5)]
    half=A(F(1,2))
    ds=directions+caps
    # All four bounds are exact support upper bounds for K0.
    h0=[g.dot(n,ta) for n in directions]+[half,half]
    posefile=comp/'finite-b-core'/'host-poses.json'
    poses=json.loads(posefile.read_bytes())['poses']
    assert len(poses)==150 and {r['index'] for r in poses}==set(range(150))
    inputs={name:json.loads((comp/name).read_bytes()) for name in {r['source'] for r in poses}}
    msource='ah-tapered-mirror-counterexample-input.json'
    base=Support(inputs[msource],arb.pi())
    mirrors={'Support':base,'Polygon':Polygon(base)}
    maxima=h0.copy()
    source_max=['K0']*4
    for r in poses:
        name=r['source']
        if 'w' in r:
            i=r['body_index']
            if name=='nonideal-pentagon-library.json':
                item=inputs[name]['proofs'][i]
                a,b=item['w'],item['deformation']
            else:
                assert name=='nonmonotone-pentagon-inputs.json'
                a,b=inputs[name][i]
            assert (F(a),F(b))==(F(r['w']),F(r['e']))
            _,support=body(ns,F(a),F(b))
        elif name==msource:
            support=mirrors[r['body_type']].support_direction
        else:
            support=Profile(inputs[name][r['body_index']]).support
        for j,n in enumerate(ds):
            height=support(old['preimage'](n,r['turn'],r['reflection']))
            if height.upper()>maxima[j].upper():
                maxima[j]=height
                source_max[j]=r['index']
            else:
                # Compare rational-point upper endpoints, including equality.
                assert height.upper()<=maxima[j].upper()
    w,e=F('.7982907'),F('.032688')
    rw,re=F('1e-10'),F('1e-11')
    roots=verify(ns,w,e,rw,re)
    _,support=body(ns,ns['box'](w-rw,w+rw),ns['box'](e-re,e+re))
    rootposes=[(0,1),(240,-1)]
    gaps=[(support(old['preimage'](n,*p))-bound.upper()).lower() for n,p,bound in zip(directions,rootposes,maxima)]
    assert min(gaps)>A(F('.000007'))
    capgaps=[]
    for k,s,i,lo in roots['cap_eliminations']:
        j=2+(i==5)
        capgaps.append((arb(lo)+half-maxima[j].upper()).lower())
    assert min(capgaps)>A(F('.0001'))
    # Absolute thickening changes every separating height, including the four non-hexagon rows.
    radius=A(F('2e-9'))
    theta=arb.pi()/18000
    assert 3*radius<arb(roots['interior_margin']).lower()
    slope=arb(roots['endpoint_margin']).lower()
    distance=(3*radius*theta/(slope*theta.sin())).upper()
    assert distance<theta
    drift=distance/arb(3).sqrt()
    motion=drift+(2/arb(3).sqrt())*(radius+drift)
    total=(radius+motion).upper()
    assert total<min(gaps+capgaps)
    files=list(inputs)+['finite-b-core/host-poses.json','deformed_reuleaux_support.py','certify_nonideal_pentagon_pose_library.py','certify_ah_tapered_minkowski_closure.py','certify_ah_tapered_mirror_counterexample.py','certify_tapered_forced_support.py']
    return dict(complete=True,python=sys.version.split()[0],flint=flint.__version__,precision=ctx.prec,
        w=str(w),e=str(e),w_radius=str(rw),e_radius=str(re),directions_degrees=angles,
        support_upper_150_plus_core=list(map(lambda h:str(h.upper()),maxima)),maximum_attainer=source_max,
        pose_order=rootposes,pose_gap_lower=list(map(str,gaps)),other_ten_gap_lower=str(min(capgaps)),
        all_angle_certificate=roots,euclidean_thickening=dict(radius='2e-9',angle_radians_upper=str(distance),
        height_change_including_host=str(total)),source_sha256={n:sha256((comp/n).read_bytes()).hexdigest() for n in sorted(files)})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('computation',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();obj=run(a.computation);a.out.write_text(json.dumps(obj,indent=2)+'\n')
    print(obj['pose_gap_lower'],obj['other_ten_gap_lower'])
