"""Tapered paired expansion, the true overlap integral, and multiple poses of triangular perturbations."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import arb,arb_series,ctx,__version__ as flint_version
from certify_ah_released_tangent import load_geometry,triangle_support


def run(archive):
    import audit_optimized_eh_chord as g
    from audit_eh_continuous_placement import polygon
    from certify_eh_upper_chain import front
    from certify_eh_joint_parameters import specialize
    ns=load_geometry(archive)
    A=ns['A']
    specialize(F('1.5497'))
    ctx.prec=160
    ctx.cap=10
    read=lambda name:json.loads((archive/name).read_bytes())
    a,b,peak,amp=map(F,('.18','1','.44','.00004'))
    p,delta=F('.950'),F('.011')
    W=g.plus(ns['z'],g.scale(A(p),ns['D']))
    d=ns['direct']
    u0=g.neg(g.rot(g.plus(g.scale(A(delta).cos(),d),g.scale(A(delta).sin(),g.rot(d)))))
    h0=g.dot(u0,W)
    def offset(w):
        s=max(0,min(1,(w-a)/(peak-a)))
        return amp*(1-w)/F('.55')*s*s*(3-2*s)
    def q(z):
        return ns['reflection'](ns['points'](z)[0][2])
    def tangent(w):
        v=q(arb_series([A(w),1]))
        x=tuple(y[0] for y in v)
        dx=tuple(y[1] for y in v)
        u=g.scale(-1/g.norm(dx),g.rot(dx))
        return x,u,g.dot(u,x)
    rows=[]
    first,*_=tangent(a)
    last,*_=tangent(b)
    for i in range(65):
        w=a+(b-a)*i/64
        _,u,h=tangent(w)
        rows.append((f'taper_{i}',u,h+A(offset(w))))
    xx=[first,*[g.line_pair(n,m,h,k) for (_,n,h),(_,m,k) in zip(rows,rows[1:])],last]
    assert all(g.cross(n,m)>0 for (_,n,_),(_,m,_) in zip(rows,rows[1:]))
    assert all(g.dot(g.rot(n),g.minus(y,x))>0 for x,y,(_,n,_) in zip(xx,xx[1:],rows))
    turns=[g.cross(g.minus(y,x),g.minus(z,y)) for x,y,z in zip(xx,xx[1:],xx[2:])]
    assert min(turns)>0
    hexn,tilt,*_=g.setup()
    pclear=min((A('.5')-g.dot(u,x)).lower() for x in xx for u in (*hexn,tilt[3],tilt[5]))
    assert pclear>0
    signs=[h0-g.dot(u0,x) for x in xx]
    keep=next(i for i,s in enumerate(signs) if s>0)
    assert all(s<0 for s in signs[:keep]) and all(s>0 for s in signs[keep:])
    label,uc,hc=rows[keep-1]
    crosspoint=g.line_pair(u0,uc,h0,hc)
    assert g.dot(g.minus(crosspoint,xx[keep-1]),g.minus(xx[keep],crosspoint))>0
    def green(l,r,origin):
        def taylor(z):
            x,y=q(arb_series([z,1]))
            return (x-origin[0])*y.derivative()-(y-origin[1])*x.derivative()
        total=A(0)
        for k in range(128):
            v0=l+(r-l)*k/128
            v1=l+(r-l)*(k+1)/128
            step=A((v1-v0)/2)
            cc=taylor(A((v0+v1)/2))
            rem=taylor(A(v0).union(A(v1)))[8]
            assert rem.is_finite()
            total+=sum(cc[i]*step**(i+1)/(i+1) for i in range(0,8,2))+rem*step**9/9
        return total
    def polygon_green(pts,origin):
        return sum(g.cross(g.minus(x,origin),g.minus(y,origin))/2 for x,y in zip(pts,pts[1:]))
    totaladded=polygon_green(xx,first)-green(a,b,first)
    augmented_cap=green(F(0),a,W)+polygon_green([*xx[:keep],crosspoint],W)
    old=read('ah-flexible-chord-certificate.json')
    ahcap=lambda j:A(F(j[0])).union(A(F(j[1])))
    oldcap=ahcap(old['shape']['cap_area_exact_interval'])
    oldarea=ahcap(read('flexible-006-eleven-eh.json')['area_exact_interval'])
    # The old-host deep cap is identical to the certified p=.95, delta=.011 calculation.
    long=json.loads(Path(__file__).with_name('ah-released-tangent-long.json').read_bytes())
    assert F(long['p'])==p and F(long['angle_radians'])==delta
    deleted=arb(long['deleted_area'])
    overlap=augmented_cap-(oldcap+deleted)
    net=deleted-totaladded+overlap
    area=oldarea-net
    assert 0<overlap<totaladded and net>A('.00001')
    # The new chord intersects the chain strictly after q(a), so the earlier paired graph and
    # this part of the chain lie entirely on the removed side. The active chord + remaining chain
    # and the retained circular/linear boundaries define the entire convex host.
    assert h0-g.dot(u0,ns['z'])<0 and h0-g.dot(u0,first)<0
    inner=polygon('0.0005522',128)[0]
    others=[x for x in inner if not x[0].startswith('chord_Z_')]
    eh=front(read('eh-two-line-cover-certificate.json'),read('eh-upper-chain-parameters.json'))
    otherEH=[(f'EH_{i}',u,h) for i,(u,h) in enumerate(zip(eh['normals'],eh['offsets']))]
    active=others+[('deep',u0,h0)]+rows[keep-1:]+otherEH
    hexnames={'B1-Q','P-F1','F1-H','E3-D1','D1-G','C3-B1'}
    for name,u,h in active:
        if name in hexnames:
            assert (h-A('.5')).contains(0)
            assert any(all((x-y).contains(0) for x,y in zip(u,z)) for z in hexn)
    sprague=('chord_Tp_','chord_G_','chord_Te_','chord_H_','chord_Q_')
    poses=[(240,-1),(0,1),(180,-1)]
    diskfree=[h-A('.5') for name,u,h in active if name not in hexnames|{'I-E3','K-C3'} and not name.startswith(sprague)]
    assert min(diskfree)>0
    def preimage(u,k,s):
        co,si=g.unit(-A(k)*arb.pi()/180)
        z=(u[0],s*u[1])
        return g.plus(g.scale(co,z),g.scale(si,g.rot(z)))
    tabs=[]
    tabs_positive=[]
    for k,s in poses:
        tab=[]
        tab_pos=[]
        for name,u,h in active:
            if name in hexnames or name.startswith(sprague) or name=='taper_64':
                continue
            pre=preimage(u,k,s)
            tri=triangle_support(g,A,ns['reflection'](pre))
            tab.append((name,pre,h,tri))
            tab_pos.append((name,pre,h,triangle_support(g,A,pre)))
        tabs.append(tab)
        tabs_positive.append(tab_pos)
    # The terminal paired tangent at q(1) is identical to the first tangent of the G arc.
    # The unit-pair map q(1) lies on that circle and the tangent join was certified in R34.
    _n,_m,_c3,_b3,*_=g.setup()
    f3=g.line_pair(_n[3],_m[2])
    gc=g.circle_line(f3,_n[5],1)
    radial=g.minus(last,gc)
    assert all((x-y).contains(0) for x,y in zip(radial,rows[-1][1]))
    assert (g.dot(radial,radial)-1).contains(0)
    assert {'I-E3','K-C3'}<=set(name for name,*_ in tabs[0])
    def cell_slacks(l,r,t0,t1,indices,js):
        mid=A((l+r)/2)
        vv=ns['points'](mid)[0]
        verts=ns['centered_points'](l,r)
        dv=ns['point_jets'](A(l).union(A(r)))[1]
        edges=[g.minus(verts[k],verts[m]) for k,m in ns['orders']]
        rad=A((r-l)/2)
        ans=[]
        for pose,ids in zip(indices,js):
            uncertain=[]
            for ix in ids:
                name,u,h,tri=tabs[pose][ix]
                cr=[g.cross(e,u) for e in edges]
                values=[(g.dot(vv[k],u)+c+rad*abs(g.dot(dv[k],u))).upper() for i,(k,c) in enumerate(ns['sector']) if not cr[i]<0 and not cr[(i+1)%10]>0]
                assert values
                val=max(values)
                if not all(h>((1-t)*val+t*tri) for t in (A(t0),A(t1))):
                    uncertain.append(ix)
            ans.append(uncertain)
        return ans
    # Verify that the exact unit-distance formulas have the same ten normal cones
    # and the same three hexagon contacts over the entire proposed w interval.
    for cell in range(512):
        w0=F('.427')+(F(1)-F('.427'))*cell/512
        w1=F('.427')+(F(1)-F('.427'))*(cell+1)/512
        v=ns['centered_points'](w0,w1)
        edges=[g.minus(v[k],v[l]) for k,l in ns['orders']]
        for k,e in enumerate(edges):
            f=edges[(k+1)%10]
            assert (e[1]>0 if k<5 else e[1]<0)
            assert g.cross(e,f)>0 and g.dot(e,f)>A('.75')
        for k in range(5):
            for l in range(k+1,5):
                if (k,l) not in ns['orders'] and (l,k) not in ns['orders']:
                    sep=g.minus(v[k],v[l])
                    assert g.dot(sep,sep)<1
        for k,u in enumerate(hexn):
            cr=[g.cross(e,u) for e in edges]
            possible={ns['sector'][j] for j in range(10) if not cr[j]<0 and not cr[(j+1)%10]>0}
            assert possible=={((1,3,0)[k%3],int(k<3))}
    # Cache discharged rows down each branch. Each successful leaf selects one rigid pose
    # for every support direction and every (w,t) in the complete leaf.
    # Preimages of the two genuine normalized cap normals are antipodes. This also
    # follows algebraically from the 60 degree pose rotation and 120 degree difference
    # in the two cap normals before reflection.
    c1=next(q for q in tabs[0] if q[0]=='K-C3')
    c2=next(q for q in tabs[2] if q[0]=='I-E3')
    assert all((x+y).contains(0) for x,y in zip(c1[1],c2[1]))
    assert (c1[2]-A('.5')).contains(0) and (c2[2]-A('.5')).contains(0)
    counts=[0]*(len(poses)+1)
    checks=0
    def sweep(l,r,tl,tr,indices,js,depth=0):
        nonlocal checks
        checks+=1
        remaining=cell_slacks(l,r,tl,tr,indices,js)
        for i,rr in zip(indices,remaining):
            if not rr:
                counts[i]+=1
                return
        if [tabs[0][k][0] for k in remaining[0]]==['K-C3'] and [tabs[2][k][0] for k in remaining[2]]==['I-E3']:
            counts[3]+=1
            return
        if depth>=36:
            raise ValueError((float(l),float(r),float(tl),float(tr),[(poses[i],[tabs[i][k][0] for k in rr][:5]) for i,rr in zip(indices,remaining)]))
        if (r-l)/F('.52')>(tr-tl)/F(1):
            m=(r+l)/2
            sweep(l,m,tl,tr,indices,remaining,depth+1)
            sweep(m,r,tl,tr,indices,remaining,depth+1)
        else:
            m=(tr+tl)/2
            sweep(l,r,tl,m,indices,remaining,depth+1)
            sweep(l,r,m,tr,indices,remaining,depth+1)
    indices=list(range(len(poses)))
    ids=[list(range(len(tab))) for tab in tabs]
    for j in range(32):
        l=F('.427')+(F(1)-F('.427'))*j/32
        r=F('.427')+(F(1)-F('.427'))*(j+1)/32
        sweep(l,r,F(0),F(1),indices,ids)
    reflected_counts, reflected_checks=counts.copy(),checks
    tabs=tabs_positive
    counts=[0]*(len(poses)+1)
    checks=0
    for j in range(32):
        l=F('.427')+(F(1)-F('.427'))*j/32
        r=F('.427')+(F(1)-F('.427'))*(j+1)/32
        sweep(l,r,F(0),F(1),indices,ids)
    fits={}
    for name in ('widened-projected-obstruction-points.json','widened-projected-obstruction-18-points.json'):
        z=read(name)
        pts=[tuple(A(F(i,z['denominator'])) for i in x) for x in z['points']]
        fits[name]=str(min((h-g.dot(u,x)).lower() for _,u,h in active for x in pts))
        assert arb(fits[name])>0
    return dict(p=str(p),delta=str(delta),start=str(a),end=str(b),peak=str(peak),amplitude=str(amp),positions=65,
                total_added=str(totaladded),old_host_deleted=str(deleted),overlap=str(overlap),surviving_added=str(totaladded-overlap),
                net_gain_R199=str(net),true_area=str(area),all_corners_P_clearance=str(pclear),min_local_cross=str(min(turns)),
                first_retained_segment=label,poses=poses,negative_triangle_leaves=reflected_counts,negative_visited=reflected_checks,positive_triangle_leaves=counts,positive_visited=checks,disk_nontrivial_slack=str(min(diskfree)),w_interval=['.427','1'],
                t_interval=['0','1'],full_Minkowski_simplex=True,historical_witnesses=fits,universal_placement_proved=False)

if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('archive',type=Path)
    cli.add_argument('output',type=Path)
    opt=cli.parse_args()
    path=opt.archive.resolve()
    sys.path.insert(0,str(path))
    obj=run(path)
    obj['script_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    obj['long_true_area_file_sha256']=sha256(Path(__file__).with_name('ah-released-tangent-long.json').read_bytes()).hexdigest()
    obj['flint']=flint_version
    obj['python']=sys.version
    opt.output.write_text(json.dumps(obj,indent=2)+'\n')
    print(obj)
