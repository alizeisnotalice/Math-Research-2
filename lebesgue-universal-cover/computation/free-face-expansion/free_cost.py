"""Free convex extensions of the single-face continuous host, with six area-implied face bounds."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys

from flint import arb, ctx, __version__ as flint_version


def run(args):
    comp, source, diff = args.comp.resolve(), args.single.resolve(), args.difference
    sys.path.insert(0,str(source))
    from certify_counterexample import geometry
    from certify_monotone_cost import quantize_chart, posed_jets, support_pins
    ns, body = geometry(comp)
    from certify_ah_tapered_minkowski_closure import geometry as prepare
    from certify_nonmonotone_ah_core import hull, double_area, SCALE, LENGTH2_LO
    from certify_tapered_forced_support import rational as A

    def read(p): return json.loads(p.read_bytes())
    def digest(p): return sha256(p.read_bytes()).hexdigest()

    old,_ = prepare(comp)
    ctx.prec=192
    g = old['g']
    ta=g.plus(old['ns']['z'],old['ns']['D'])
    v=g.minus(old['q'](arb(1)),ta)
    chart,outward,inward=quantize_chart(ns,g,ta,v,SCALE)
    x0=g.neg(g.rot(v))
    assert g.dot(v,v)>A(LENGTH2_LO)
    normals=[]
    for j in range(-args.steps,args.steps+1):
        p=g.plus(x0,g.scale(arb(j)/(3*args.steps),v))
        normals.append(g.scale(1/g.norm(p),p))

    host = source/'single-H1-input.json'
    base_file = comp/'single-repair'/'single-repair.json'
    oldbase=F(read(base_file)['base_area'][0])
    d=read(diff)
    hostbounds=list(map(F,[d['single']['lower'],d['single']['upper']]))
    assert digest(host)==d['input_sha256']['single-H1-input.json']
    cap_file=source/'cost-86.json'
    cert=read(cap_file)
    assert cert['provenance']['host']==digest(host)
    assert cert['provenance']['script']==digest(source/'certify_monotone_cost.py')
    cap_upper=F(cert['original_AH_cap_area_upper'])
    gain=args.gain
    assert 0<gain
    threshold=(cap_upper+gain)*2*SCALE**2/LENGTH2_LO

    n,m,c3,b3,*_ = g.setup()
    half=arb(1)/2
    def project(x,n): return g.plus(x,g.scale(half-g.dot(n,x),n))
    b1=g.line_pair(n[0],n[1]);f1=g.line_pair(n[2],n[3]);d1=g.line_pair(n[4],n[5])
    e3=g.line_pair(n[4],m[3]);f3=g.line_pair(n[3],m[2]);G=g.circle_line(f3,n[5],1)
    faces=[(c3,b1),(b1,project(e3,n[1])),(project(G,n[2]),f1),(f1,project(c3,n[3])),(e3,d1),(d1,G)]
    for i,seg in enumerate(faces):
        assert all((g.dot(n[i],z)-half).contains(0) for z in seg)
        assert all(h>g.dot(t,z) for _,t,h in old['otherEH'] for z in seg)
        assert all(g.dot(x0,g.minus(z,ta))<0 for z in seg)
    lengths=[g.norm(g.minus(a,b)) for a,b in faces]
    u=[2*A(gain if j!=1 else gain+hostbounds[1]-oldbase)/lengths[j] for j in range(6)]
    assert max(z.upper() for z in u)<arb('.1')
    plus_length=min(lengths[j].lower() for j in (1,3,5))
    minus_length=min(lengths[j].lower() for j in (0,2,4))
    q_upper=2*A(gain+hostbounds[1]-oldbase)/plus_length
    q_lower=-2*A(gain)/minus_length
    print('lengths',list(map(float,lengths)),'q',float(q_lower),float(q_upper),flush=True)
    dx= (-u[3]).union(u[0]); dz=(-u[5]).union(u[2])
    rt=arb(3).sqrt()
    capnorm=[old['tilt'][j] for j in (3,5)]
    cappairs=[(project(b3,m[3]),e3),(project(f3,m[5]),c3)]
    for seg in cappairs:
        assert g.norm(g.minus(*seg))>A(F('.2535'))
        assert all(h>g.dot(t,z) for _,t,h in old['otherEH'] for z in seg)
        assert all(g.dot(x0,g.minus(z,ta))<0 for z in seg)
    points=[(0,0),(SCALE,0)]
    for name in ('pentagons.json','curves.json'):
        tab=read(comp/'finite-b-core'/name)
        vertices=tab.get('combined_hull') or tab['vertices_ccw_integer']
        den=tab.get('scale') or tab['grid']
        assert F(tab.get('coordinate_rounding_error_strict') or tab['coordinate_error_asserted'])<=F(2,den)
        err=arb(0,arb(2)/den)
        for x,y in vertices:
            p=inward((arb(x)/den+err,arb(y)/den+err))
            if p is not None: points.append(p)
    inner=hull(points)
    witness=Path(args.witness).resolve() if args.witness else source/"cost-input.json"
    prob=read(witness)
    vs,_=body(ns,F(prob['w']),F(prob['e']))
    es=[ns['sub'](vs[i],vs[j]) for i,j in ns['orders']]
    trig=ns['trig']; roots=[trig(ns['A'](k)) for k in (0,60,120)]
    counts={}; minimum=None
    for s in (1,-1):
        stack=[(F(i),F(i+1),0) for i in range(360)]; stat={}; max_depth=0
        while stack:
            lo,hi,k=stack.pop()
            mid=(lo+hi)/2; omega=ns['box'](lo,hi)
            hw=arb(0,ns['A']((hi-lo)/2)*ns['angle_unit'])
            fn,tc,q0=support_pins(ns,vs,es,s,ns['A'](mid))
            trig_interval=trig(-omega)
            def deriv(n):return -ns['supports'](vs,es,trig_interval,(n[0],s*n[1]),True)
            d0,d1,d2=map(deriv,roots)
            q=q0+hw*(d1-d0-d2)
            vel=(-d0,(-d0-2*d2)/rt)
            def close(label):stat[label]=stat.get(label,0)+1
            if q<q_lower or q>q_upper or q+dx+dz < -u[4] or q+dx+dz > u[1]:
                close('area_implied_hexagon');continue
            # Propagate all three width equations into the two free translations.
            bounds=(-u[4]).union(u[1])-q
            x=dx.intersection(bounds-dz);z=dz.intersection(bounds-dx)
            residual=(x,(x+2*z)/rt)
            if any((fn(n)+g.dot(n,tc)-half+hw*(deriv(n)+g.dot(n,vel))+g.dot(n,residual))*A(F('.2535'))/2>A(gain) for n in capnorm):
                close('slanted_face_area');continue
            if hi-lo<=F(1,1000):
                pinned=g.plus(tc,g.scale(hw,vel))
                def test_cell(x,z):
                    nonlocal minimum
                    residual=(x,(x+2*z)/rt)
                    off=g.plus(pinned,residual)
                    material=[inward(posed_jets(ns,vs,es,omega,s,norm,off)) for norm in normals]
                    area2=double_area(hull(inner+tuple(p for p in material if p is not None)))
                    ah=F(area2,2*SCALE**2)*LENGTH2_LO
                    if area2>threshold:
                        minimum=area2 if minimum is None else min(minimum,area2)
                        return 'AH_area'
                    height=(x,q+x+z,z,-x,-q-x-z,-z)
                    credits=[]
                    for j,delta in enumerate(height):
                        value=arb(0)
                        if delta>0:
                            point=posed_jets(ns,vs,es,omega,s,n[j],off)
                            if g.dot(x0,g.minus(point,ta))<0:
                                value=delta*lengths[j]/2
                        credits.append(value.lower())
                    even=sum(credits[::2]);odd=sum(credits[1::2])
                    credit=max(even.lower(),(odd-A(hostbounds[1]-oldbase)).lower(),(credits[3]+credits[5]).lower())
                    if A(ah-cap_upper)+credit>A(gain): return 'AH_plus_faces'
                    return None
                kind=test_cell(x,z)
                if kind:
                    close(kind);continue
                if hi-lo<=F(1,65536):
                    def split(y):
                        return (y.lower().union(y.mid()),y.mid().union(y.upper()))
                    def certify_translation(x,z,depth):
                        if q+x+z < -u[4] or q+x+z > u[1]:return True
                        feasible=(-u[4]).union(u[1])-q
                        x=x.intersection(feasible-z);z=z.intersection(feasible-x)
                        if test_cell(x,z):return True
                        if depth>=args.free_depth:return False
                        if float(x.rad())>=float(z.rad()):
                            return all(certify_translation(xx,z,depth+1) for xx in split(x))
                        return all(certify_translation(x,zz,depth+1) for zz in split(z))
                    if certify_translation(x,z,0):
                        close('joint_translation_partition');continue
            if k>=22:
                raise ValueError((s,str(lo),str(hi),str(q),stat))
            stack.extend(((lo,mid,k+1),(mid,hi,k+1)));max_depth=max(max_depth,k+1)
        counts[str(s)]=dict(leaves=stat,max_depth=max_depth)
        print(s,counts[str(s)],flush=True)
    return dict(strict_gain=str(gain),strict_area_lower=str(hostbounds[0]+gain),witness=prob,
                old_supporting_segment_lengths=[str(z) for z in lengths],displacement_upper=[str(z) for z in u],
                cap_upper=str(cap_upper),minimum_closed_AH_area=str(F(minimum,2*SCALE**2)*LENGTH2_LO) if minimum else None,
                AH_normals=len(normals),angle_cover=counts,python=platform.python_version(),flint=flint_version,precision=ctx.prec,
                material_sources={str(p.name):digest(p) for p in (source/'certify_counterexample.py',source/'certify_monotone_cost.py',witness,host,cap_file,base_file,diff,Path(__file__))})

if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('comp',type=Path);p.add_argument('single',type=Path);p.add_argument('difference',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--gain',type=F,default=F('3e-8'));p.add_argument('--steps',type=int,default=160)
    p.add_argument('--witness',type=Path);p.add_argument('--free-depth',type=int,default=6)
    args=p.parse_args()
    args.output.write_text(json.dumps(run(args),indent=2)+'\n')
