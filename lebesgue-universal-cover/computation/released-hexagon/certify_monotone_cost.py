"""Lower-bound the area required by any H1-released convex enlargement of C153."""
import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import platform
import sys
import time

from flint import arb, ctx, __version__ as flint_version
sys.path.insert(0, str(Path(__file__).resolve().parent))
from certify_counterexample import digest, geometry


def quantize_chart(ns, g, ta, v, q):
    v2=g.dot(v,v)
    def chart(x):
        d=g.minus(x,ta)
        return g.dot(d,v)/v2,-g.cross(v,d)/v2
    def outward(p):
        x,y=chart(p)
        i=int(float(x.mid())*q)
        j=int(float(y.mid())*q)
        lowx=int(float((x*q).lower()))-4; hix=int(float((x*q).upper()))+4
        lowy=int(float((y*q).lower()))-4; hiy=int(float((y*q).upper()))+4
        assert x*q>lowx and x*q<hix and y*q>lowy and y*q<hiy
        return tuple(product((lowx,hix),(lowy,hiy)))
    def inward(p):
        x,y=chart(p)
        if not y>arb(100)/q:
            return None
        # The common cross section is [Y * sup(x/y), 1 + Y * inf((x-1)/y)].
        # Pick the highest grid level that leaves a rational point between them.
        left=(x/y).upper()
        right=((x-1)/y).lower()
        width_loss=left-right
        ymax=y.lower()
        if width_loss>0:
            ymax=min(ymax,(1/width_loss).lower())
        iy=int(float(ymax)*q)-16
        if iy<=0:return None
        for skip in (0,32,256):
            Y=arb(iy-skip)/q
            if not 0<Y<y:continue
            l,r=Y*left,1+Y*right
            m=int(float(((l+r)*q/2).mid()))
            if l*q<m<r*q:
                return m,iy-skip
        return None
    return chart,outward,inward


def posed_jets(ns, vs, es, alpha, s, normal, offset):
    trig=ns['trig']; dot=ns['dot']; sub=ns['sub']; add=ns['add']; mul=ns['mul']; cross=ns['cross']
    co,si=trig(alpha)
    nx,ny=normal
    u=(co*nx+s*si*ny,-si*nx+s*co*ny)
    votes=[cross(e,u) for e in es]
    candidates=[]
    for i,k in enumerate(votes):
        if k<0 or votes[(i+1)%len(es)]>0:continue
        j,r=ns['sector'][i]
        p=vs[j]
        rotated=(co*p[0]-si*p[1],s*(si*p[0]+co*p[1]))
        candidates.append(add(add(rotated,mul(r,normal)),offset))
    assert candidates
    x,y=candidates[0]
    for p,q in candidates[1:]:
        x=x.union(p);y=y.union(q)
    return x,y


def support_pins(ns, vs, es, s, alpha):
    A=ns['A']; trig=ns['trig']
    half=arb('.5');rt3=arb(3).sqrt()
    co,si=trig(-alpha)
    def raw(n):return ns['supports'](vs,es,(co,si),(n[0],s*n[1]))
    h0,h1,h2=[raw(trig(A(i))) for i in (0,60,120)]
    a,b=half-h0,half-h2
    return raw,(a,(2*b+a)/rt3),half+h1-h0-h2


def clip_y(poly,y):
    ans=[]
    for a,b in zip(poly[-1:]+poly[:-1],poly):
        sa=a[1]>=y;sb=b[1]>=y
        if sa!=sb:
            t=(y-a[1])/(b[1]-a[1]);ans.append((a[0]+t*(b[0]-a[0]),y))
        if sb:ans.append(b)
    return ans


def main(args):
    start=time.monotonic();comp=args.comp.resolve();ns,body=geometry(comp)
    from certify_ah_tapered_minkowski_closure import geometry as old_geometry
    from certify_nonmonotone_ah_core import hull,double_area, BASE_LO,LENGTH2_LO,SCALE
    from certify_tapered_forced_support import rational
    from flint import arb_mat
    old,_=old_geometry(comp);g=old['g'];ctx.prec=192
    ta=g.plus(old['ns']['z'],old['ns']['D']);tp=old['q'](arb(1));v=g.minus(tp,ta);v2=g.dot(v,v)
    assert v2>rational(LENGTH2_LO)
    chart,outward,inward=quantize_chart(ns,g,ta,v,SCALE)
    cfg=json.loads(args.host.read_bytes()); assert cfg['translation_rule']=='pin_hexagon_faces_0_and_2'
    core,pent,curves=[json.loads((comp/'finite-b-core'/n).read_bytes()) for n in ('core-outer.json','pentagons.json','curves.json')]
    inner=[(0,0),(SCALE,0)]; outer=[]
    for x,y in core['outer_integer_vertices']:
        outer.extend(outward((arb(x)/core['vertex_grid'],arb(y)/core['vertex_grid'])))
    assert pent['grid']==curves['scale'] and pent['directions']==curves['N']
    old_points=pent['vertices_ccw_integer']+curves['combined_hull']; Den=pent['grid'];eps=arb(2)/Den
    assert F(pent['coordinate_rounding_error_strict'])<=F(2,Den)
    assert F(curves['coordinate_error_asserted'])<=F(2,Den)
    for x,y in old_points:
        p=(arb(x)/Den,arb(y)/Den)
        outer.extend(outward(p))
        c=inward(tuple(z+arb(0,eps) for z in p))
        if c is not None:inner.append(c)
    N=pent['directions'];planes=[ns['trig'](ns['A'](F(360*j,N)))for j in range(N)]
    for w,e,z,s in cfg['bodies']:
        vs,_=body(ns,F(w),F(e)); es=[ns['sub'](vs[i],vs[j]) for i,j in ns['orders']]
        a=ns['A'](F(z));_,t,_=support_pins(ns,vs,es,s,a)
        for n in planes:
            p=posed_jets(ns,vs,es,a,s,n,t)
            outer.extend(outward(p))
            c=inward(p)
            if c is not None:inner.append(c)
    li,ri=cfg['family_interval']['first'],cfg['family_interval']['last']
    l,r=cfg['bodies'][li],cfg['bodies'][ri]
    assert l[0]==r[0] and l[2:]==r[2:] and F(l[1])<F(r[1])
    interpolation=ns['A'](12)*ns['A'](F(r[1])-F(l[1]))**2/8
    delta=(2*(arb.pi()/(2*N)).sin()**2+arb(3)/Den+interpolation)/v2.sqrt()
    R=F('.000000047'); assert delta<rational(R)
    poly=hull(outer)
    clip=clip_y([(F(x,SCALE),F(y,SCALE)) for x,y in poly],-R)
    chart_area=abs(double_area(tuple(clip)))/2
    perimeter=sum(((rational(x-u)**2+rational(y-vv)**2).sqrt() for (x,y),(u,vv) in zip(clip,clip[1:]+clip[:1])),arb(0))
    excess_upper=(rational(chart_area)+rational(R)*perimeter+arb.pi()*rational(R)**2)*v2
    inner=hull(inner)
    inner_a2=double_area(inner)
    inner_area=F(inner_a2,2*SCALE**2)*LENGTH2_LO
    print('C153 cap area lower',float(inner_area),'upper',float(excess_upper.upper()),'time',round(time.monotonic()-start,2),flush=True)
    relative=F(args.gain)
    qlo=F(int(float(excess_upper.upper())*10**16)+10,10**16)
    assert rational(qlo)>excess_upper
    threshold=(qlo+relative)*2*SCALE**2/LENGTH2_LO
    inp=json.loads(args.test.read_bytes());vs,_=body(ns,F(inp['w']),F(inp['e']));es=[ns['sub'](vs[i],vs[j]) for i,j in ns['orders']]
    assert all(ns['dot'](p,p)<arb(1)/3 for p in vs)
    caps=[old['tilt'][j] for j in (3,5)]
    roots=[ns['trig'](ns['A'](i))for i in (0,60,120)]
    # The original core's straight cap sides remain present after both cuts.
    n,m,c3,b3,*_=g.setup()
    def project(p,normal):
        return g.plus(p,g.scale(arb(1)/2-g.dot(p,normal),normal))
    e3=g.line_pair(n[4],m[3]); f3=g.line_pair(n[3],m[2])
    segments=[(project(b3,m[3]),e3),(project(f3,m[5]),c3)]
    DD=F('0.2535')
    for p,q in segments:
        assert g.norm(g.minus(p,q))>rational(DD)
        assert all(h>g.dot(t,z) for _,t,h in old['otherEH'] for z in (p,q))
        assert all(g.dot(g.neg(g.rot(v)),g.minus(z,ta))<0 for z in (p,q))
    old_caps=json.loads((comp/'finite-b-core'/'obstruction.json').read_bytes())
    assert old_caps['complete'] and all(arb(h)<=arb('.5') for h in old_caps['support_upper_150_plus_core'][2:])
    for w,e,z,s0 in cfg['bodies']:
        vs0,_=body(ns,F(w),F(e));es0=[ns['sub'](vs0[i],vs0[j]) for i,j in ns['orders']]
        fn,t0,q0=support_pins(ns,vs0,es0,s0,ns['A'](F(z)))
        assert all(fn(u)+g.dot(u,t0)+interpolation<arb('.5') for u in caps)
    supports=ns['supports']
    trig=ns['trig']
    normals=[];x0=g.neg(g.rot(v))
    for j in range(-args.steps,args.steps+1):
        p=g.plus(x0,g.scale(arb(j)/(3*args.steps),v)); normals.append(g.scale(1/g.norm(p),p))
    vals={}; best_area=None
    for s in (1,-1):
        stack=[(F(i),F(i+1),0)for i in range(360)];counts={};depth=0
        while stack:
            l,r,k=stack.pop();mid=(l+r)/2
            sig=ns['box'](l,r); rads=ns['A']((r-l)/2)*ns['angle_unit']
            fn,t0,q0=support_pins(ns,vs,es,s,ns['A'](mid))
            csa=trig(-sig)
            def derivative(u):return -supports(vs,es,csa,(u[0],s*u[1]),True)
            d0,d1,d2=map(derivative,roots)
            interval=arb(0,rads)
            q=q0+interval*(d1-d0-d2)
            dtx=-d0; dty=(-2*d2-d0)/arb(3).sqrt()
            t=(t0[0]+interval*dtx,t0[1]+interval*dty)
            row=None
            if q<0 or q>rational(args.eta):row='strip'
            if row is None:
                limits=[fn(n)+ns['dot'](n,t0)-arb('.5')+interval*(derivative(n)+ns['dot'](n,(dtx,dty))) for n in caps]
                if any(d*rational(DD)/2>rational(relative) for d in limits):row='cap'
            if row is None and r-l<=F(1,1000):
                choices=[]
                for n in normals:
                    z=posed_jets(ns,vs,es,sig,s,n,t)
                    p=inward(z)
                    if p is not None:choices.append(p)
                aug=hull(inner+tuple(choices));area2=double_area(aug)
                if area2*threshold.denominator>threshold.numerator:
                    row='AH'; lo=F(area2,2*SCALE**2)*LENGTH2_LO
                    best_area=lo if best_area is None else min(lo,best_area)
            if row is None:
                assert k<22,(s,str(l),str(r),str(q),str(relative),counts)
                stack.extend(((l,mid,k+1),(mid,r,k+1)))
                continue
            counts[row]=counts.get(row,0)+1;depth=max(depth,k)
        print('proved',s,counts,'depth',depth,'seconds',round(time.monotonic()-start,2),flush=True)
        vals[str(s)]=dict(counts=counts,max_depth=depth)
    result=dict(test=inp,eta=str(args.eta),strict_relative_gain=str(relative),original_AH_cap_area_lower=str(inner_area),
                original_AH_cap_area_upper=str(qlo),new_AH_alternatives_cap_area_lower=str(best_area),
                original_core_area_lower=str(BASE_LO),absolute_area_baseline_plus_relative_gain=str(BASE_LO+inner_area+relative),
                chart_grid=SCALE,area_target_squared_chord_lower=str(LENGTH2_LO),core_inner_polygon_vertices=len(inner),outer_clipped_vertices=len(clip),
                profile_release=vals,genuinely_pinned_original_cap_length=str(DD),chart_radius_upper=str(R),normal_count=len(normals),
                provenance=dict(host=digest(args.host),candidate=digest(args.test),script=digest(Path(__file__)),
                    python=platform.python_version(),flint=flint_version,precision=ctx.prec))
    args.out.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    if not __debug__:
        raise RuntimeError("Assertions must be enabled.")
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('comp',type=Path);p.add_argument('host',type=Path);p.add_argument('test',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--gain',type=F,default=F('1.8e-8'));p.add_argument('--eta',type=F,default=F('0.000000301'));p.add_argument('--steps',type=int,default=180)
    main(p.parse_args())
