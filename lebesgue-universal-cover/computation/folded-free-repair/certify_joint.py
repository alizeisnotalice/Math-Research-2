"""Check a simultaneous area lower bound over two independent all-angle poses."""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import platform
import sys
import time

from flint import arb, arb_series, ctx, __version__ as flint_version
from certify_counterexample import digest, geometry
from certify_monotone_cost import posed_jets, quantize_chart, support_pins

def dependency_log(roots):
    """Return an inventory of source and data files read under the named roots."""
    from importlib.util import source_from_cache

    recorded = set()
    def audit(event, args):
        if event != 'open' or not isinstance(args[0], (str, bytes)):
            return
        name = args[0].decode() if isinstance(args[0], bytes) else args[0]
        name = Path(name).resolve()
        if name.suffix == '.pyc':
            try:
                name = Path(source_from_cache(str(name))).resolve()
            except ValueError:
                return
        for label, root in roots.items():
            if name.is_relative_to(root) and name.is_file():
                recorded.add((label, name.relative_to(root), name))

    sys.addaudithook(audit)
    return lambda: {str(Path(label) / part): digest(path)
                    for label, part, path in sorted(recorded)}

def cap_bounds(comp,repair,fixture,cfg,ns,body,anchors,old,core):
    from certify_nonmonotone_ah_core import hull,double_area,SCALE,LENGTH2_LO
    from certify_tapered_forced_support import Profile,rational as A
    from certify_ah_tapered_mirror_counterexample import Support,Polygon

    g=old['g']
    ta,v=anchors
    cross=g.neg(g.rot(v))
    v2=g.dot(v,v)
    chart,_,inward=quantize_chart(ns,g,ta,v,SCALE)
    folder=comp/'finite-b-core'
    poses,pent,curves=[json.loads((folder/q).read_bytes()) for q in ('host-poses.json','pentagons.json','curves.json')]
    D=pent['grid']
    N=pent['directions']
    assert N==curves['N'] and D==curves['scale']
    ex=arb(2)/D
    assert F(pent['coordinate_rounding_error_strict'])<=F(2,D)
    assert F(curves['coordinate_error_asserted'])<=F(2,D)
    sampled=hull(tuple(map(tuple,pent['vertices_ccw_integer']+curves['combined_hull'])))
    pin=[]
    points=[(0,0),(SCALE,0)]
    for x,y in sampled:
        item=inward((arb(x)/D+arb(0,ex),arb(y)/D+arb(0,ex)))
        if item:points.append(item)
    wind=[g.plus(cross,g.scale(arb(i)/(3*1200),v))for i in range(-1200,1201)]
    wind=[g.scale(1/g.norm(z),z)for z in wind]
    for w,e,aa,s in cfg['bodies']:
        b,_=body(ns,F(w),F(e))
        es=[ns['sub'](b[i],b[j])for i,j in ns['orders']]
        a=ns['A'](F(aa))
        fn,t,_=support_pins(ns,b,es,s,a)
        pin.append((fn,t))
        for n in wind:
            item=inward(posed_jets(ns,b,es,a,s,n,t))
            if item:points.append(item)
    left,right=[cfg['family_interval'][x] for x in ('first','last')]
    assert (left,right)==(1,2) and len(cfg['bodies'])==3
    l,r=[cfg['bodies'][i]for i in (left,right)]
    assert l[0]==r[0] and l[2:]==r[2:]
    w,lo,hi=F(l[0]),F(l[1]),F(r[1])
    assert lo<hi
    pbox=ns['box'](lo,hi)
    body(ns,w,pbox)
    ell,t0=ns['chart'](ns['add'](ns['z'],ns['mul'](A(w),ns['D'])))
    partner=ns['psi'](ell)
    t=t0-arb_series([pbox,1])
    for x,y in (ns['X'](ell,t),ns['X'](partner,-t)):
        assert 2*(x[2]**2+y[2]**2).sqrt()<4
    # Each fixed-direction support derivative is 4-Lipschitz. Pinned translation
    # adds at most 8, since the two-row strips' inverse has norm sqrt(2).
    interpol=A(F(12,8))*A(hi-lo)**2
    assert interpol<arb('8.86e-12')
    datasets={s:json.loads((comp/s).read_bytes()) for s in {row['source']for row in poses['poses']}}
    mirror=Support(datasets['ah-tapered-mirror-counterexample-input.json'],arb.pi())
    mirrors={'Support':mirror,'Polygon':Polygon(mirror)}
    funcs=[]
    for pos in poses['poses']:
        if 'w' in pos:
            _,fn=body(ns,F(pos['w']),F(pos['e']))
        elif 'body_type' in pos:
            fn=mirrors[pos['body_type']].support_direction
        else:
            fn=Profile(datasets[pos['source']][pos['body_index']]).support
        co,si=ns['trig'](ns['A'](F(pos['turn'])))
        funcs.append((co,si,pos['reflection'],fn))
    def fan_support(x,y,co,si,s,h):
        return h((co*x+s*si*y,-si*x+s*co*y))

    Q=repair.Q
    samples=fixture['outer_boundary']
    normals=set()
    for (x,y),(u,vv) in zip(samples,samples[1:]+samples[:1]):
        d=((vv-y)**2+(x-u)**2)**.5
        if d>1.e-15:
            normals.add((round(10**12*(vv-y)/d),round(10**12*(x-u)/d)))
    assert v2>A(LENGTH2_LO)
    assert A(F(32))*v2.sqrt()>arb(3) and g.norm(ta)<arb(1)
    poly=[(-32,0),(32,0),(32,32),(-32,32)]
    rad=2*(arb.pi()/(2*N)).sin()**2+arb(3)/D
    f=10**16
    all_exact=0
    for a,b in sorted(normals):
        n=arb(a*a+b*b).sqrt()
        U=(arb(a)/n,arb(b)/n)
        core_hi=arb(max(a*x+b*y for x,y in core))/(Q*n)
        sample_hi=arb(max(a*x+b*y for x,y in sampled))/(D*n)+rad
        if sample_hi>core_hi:
            try:
                exact_hi=max(fan_support(*U,*row).upper()for row in funcs)
                if exact_hi<sample_hi:
                    sample_hi=exact_hi
                    all_exact+=1
            except AssertionError:
                pass
        old_h=max(core_hi.upper(),sample_hi.upper())
        fs=[(h(U)+g.dot(U,t)).upper()for h,t in pin]
        H=max(old_h,fs[0],(max(fs[1:])+interpol).upper())-g.dot(U,ta)
        w1,w2=g.dot(U,v),g.dot(U,cross)
        u,z=round(float((f*w1).mid())),round(float((f*w2).mid()))
        hi=H*f+32*(abs(w1*f-u)+abs(w2*f-z))
        c=round(float(hi.upper()))+8
        assert arb(c)>hi
        poly=repair.clip(poly,u,z,c)
    area=A(F(double_area(poly),2))*v2
    scale=10**18
    k=round(float((area*scale).upper()))+16
    assert arb(k)>area*scale
    return hull(points),F(k,scale),dict(inner_pts=len(hull(points)),upper_pts=len(poly),
                                            exact_refinements=all_exact,proposed_supports=len(normals))

def run(a):
    start=time.monotonic()
    comp=a.comp.resolve()
    manifest=dependency_log({'computation':comp,'single-repair':a.single_repair.resolve()})
    sys.path.insert(0,str(comp))
    sys.path.insert(0,str(a.single_repair.resolve()))
    import certify_single_repair as repair
    from certify_nonmonotone_ah_core import hull,double_area,SCALE,LENGTH2_LO
    from certify_tapered_forced_support import rational
    ns,body=geometry(comp)
    old,core,_=repair.refined_core(comp,2048,hull)
    ctx.prec=192
    g=old['g']; ta=g.plus(old['ns']['z'],old['ns']['D']); v=g.minus(old['q'](arb(1)),ta)
    chart,_,inward=quantize_chart(ns,g,ta,v,SCALE)
    cfg=json.loads(a.host.read_bytes()); assert cfg['translation_rule']=='pin_hexagon_faces_0_and_2'
    fixture=json.loads(a.fixture.read_bytes())
    test=json.loads(a.pair.read_bytes())
    gain=F(test['relative_gain']); eta=rational(F(test['eta']))
    inner,upper,quality=cap_bounds(comp,repair,fixture,cfg,ns,body,(ta,v),old,core)
    print('cap bounds',str(F(double_area(inner),2*SCALE*SCALE)*LENGTH2_LO),float(upper),quality,flush=True)
    min_area=(upper+gain)*2*SCALE*SCALE/LENGTH2_LO
    cuts=[old['tilt'][i]for i in (3,5)]
    DD=F('.2535')
    n,m,c3,b3,*_=g.setup()
    project=lambda p,n:g.plus(p,g.scale(arb(1)/2-g.dot(p,n),n))
    segs=[(project(b3,m[3]),g.line_pair(n[4],m[3])),(project(g.line_pair(n[3],m[2]),m[5]),c3)]
    for aa,bb in segs:
        assert g.norm(g.minus(aa,bb))>rational(DD)
        assert all(h>g.dot(n,p)for _,n,h in old['otherEH'] for p in (aa,bb))
    receipt=json.loads((comp/'finite-b-core'/'obstruction.json').read_bytes()); assert receipt['complete']
    assert all(arb(q)<=arb('.5')for q in receipt['support_upper_150_plus_core'][2:4])
    for w,e,k,s in cfg['bodies']:
        vs,_=body(ns,F(w),F(e))
        ee=[ns['sub'](vs[i],vs[j])for i,j in ns['orders']]
        fn,t,_=support_pins(ns,vs,ee,s,ns['A'](F(k)))
        assert all(fn(u)+g.dot(u,t)+arb('8.86e-12')<arb('.5')for u in cuts)
    bodies=[]
    for w,e,*_ in test['test_bodies']:
        xy,_=body(ns,F(w),F(e))
        bodies.append((xy,[ns['sub'](xy[i],xy[j])for i,j in ns['orders']]))
    trig=ns['trig']; rt3=arb(3).sqrt()
    base_supports=[trig(ns['A'](k))for k in (0,60,120)]
    x0=g.neg(g.rot(v)); dirs=[]
    for j in range(-a.steps,a.steps+1):
        p=g.plus(x0,g.scale(arb(j)/(3*a.steps),v));dirs.append(g.scale(1/g.norm(p),p))
    for n0,n1 in zip(fixture['outer_boundary'],fixture['outer_boundary'][1:]+fixture['outer_boundary'][:1]):
        direction=tuple(arb(str(b-a))for a,b in zip(n0,n1))
        dd=g.neg(g.rot(direction))
        if g.dot(dd,x0)>0:dirs.append(g.scale(1/g.norm(dd),dd))

    @lru_cache(maxsize=25000)
    def check(item,s,lo,hi,with_pts):
        vs,ee=bodies[item]
        mid=(lo+hi)/2
        sig=ns['box'](lo,hi)
        fn,t,q0=support_pins(ns,vs,ee,s,ns['A'](mid))
        rr=arb(0,ns['A']((hi-lo)/2)*ns['angle_unit'])
        tr=trig(-sig)
        def deriv(n):return -ns['supports'](vs,ee,tr,(n[0],s*n[1]),True)
        d0,d1,d2=map(deriv,base_supports)
        q=q0+rr*(d1-d0-d2)
        if q<0 or q>eta:return 'strip',()
        dt=(-d0,(-2*d2-d0)/rt3)
        for n in cuts:
            loss=fn(n)+g.dot(n,t)-arb('.5')+rr*(deriv(n)+g.dot(n,dt))
            if loss*rational(DD)/2>rational(gain):return 'cap',()
        if not with_pts:return None,()
        placed=(t[0]+rr*dt[0],t[1]+rr*dt[1])
        new=[p for n in dirs if (p:=inward(posed_jets(ns,vs,ee,sig,s,n,placed)))]
        p=hull(inner+tuple(new))
        return ('singleAH' if double_area(p)>min_area else None),p

    def split(t):
        s,l,r=t
        m=(l+r)/2
        return (s,l,m),(s,m,r)

    def init(i):
        allowed=(test['first_handedness'],) if i==0 and 'first_handedness' in test else (1,-1)
        stack=[(s,F(j),F(j+1))for s in allowed for j in range(360)]
        good=[];hist=Counter()
        while stack:
            x=stack.pop()
            ans,_=check(i,*x,False)
            if ans:hist[ans]+=1
            elif x[2]-x[1]<=F(1,256):good.append(x)
            else:stack.extend(split(x))
        return good,dict(hist)
    prepared=[init(k)for k in range(len(bodies))]
    panels=[r[0]for r in prepared];pruning=[r[1]for r in prepared]
    print('init',list(map(len,panels)),pruning,'seconds',round(time.monotonic()-start,2),flush=True)
    stack=list(product(*panels))
    stat=Counter(); budget=0; open_boxes=[]; proved_area=None
    while stack:
        states=stack.pop()
        budget+=1
        material=[]
        for i,sub in enumerate(states):
            answer,points=check(i,*sub,True)
            if answer:break
            material+=points
        if answer:
            stat[answer]+=1;continue
        area=double_area(hull(material))
        if area>min_area:
            stat['jointAH']+=1
            proved_area=area if proved_area is None else min(proved_area,area)
            continue
        widths=[sub[2]-sub[1]for sub in states]
        index=max(range(len(widths)),key=widths.__getitem__)
        if widths[index]>F(1,2**a.depth):
            stack.extend(states[:index]+(part,)+states[index+1:]for part in split(states[index]))
        else:
            open_boxes.append(states)
            if len(open_boxes)>=10:break
        if budget>=a.budget:break
    success=not stack and not open_boxes
    print('proved',success,stat,'open',len(open_boxes),'pending',len(stack),'secs',time.monotonic()-start,flush=True)
    result=dict(complete=success,relative_gain=str(gain),state_counts=dict(stat),starting_phase_counts=list(map(len,panels)),first_handedness=test.get('first_handedness'),
                first_pruning=pruning,cap_quality=quality,cap_area_upper=str(upper),interval_state_budget=budget,
                pending_states=len(stack),unresolved=[[list(map(str,z))for z in row]for row in open_boxes],
                all_repository_inputs=manifest(),
                deps={p.name:digest(p)for p in (Path(__file__),a.host,a.pair,a.fixture,a.single_repair/'certify_single_repair.py',
                      comp/'released-hexagon'/'certify_monotone_cost.py',comp/'released-hexagon'/'certify_counterexample.py')},
                depth=a.depth,steps=a.steps,python=platform.python_version(),flint=flint_version,precision=ctx.prec)
    if proved_area is not None:
        result['minimum_used_joint_area']=str(F(proved_area,2*SCALE*SCALE)*LENGTH2_LO)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for field in ('comp','host','pair','fixture','single_repair','output'):p.add_argument(field,type=Path)
    p.add_argument('--steps',type=int,default=300)
    p.add_argument('--depth',type=int,default=24)
    p.add_argument('--budget',type=int,default=60000)
    run(p.parse_args())
