"""Quantify compulsory penetration of the two AH supporting lines in the larger P host."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import inspect,json
from pathlib import Path
import sys
from flint import arb,arb_series,ctx

def run(comp,joint,strong):
    sys.path.insert(0,str(comp))
    from certify_sigma_paired_rigidity import rig
    from certify_low_sigma_ah_budget import g,b,model,sub,add,dot,mul,rot
    from certify_three_sigma_area_budget import exact,pack
    ctx.prec=192
    old=json.loads(joint.read_bytes());area=json.loads(strong.read_bytes())
    lam=F(1,20000);ll=exact(lam);w=F('.435')
    src=inspect.getsource(rig)
    assert src.count('if gap>0:')==1
    ns=dict(rig.__globals__);key="if gap>g.exact(F(13,10**6)):"
    exec(compile(src.replace('if gap>0:',key),'rig_positive_margin','exec'),ns)
    newrig=ns['rig']
    sig=[];harm=[];tri=[];cap=[];tot=0
    ws=[w]+[F(i,20) for i in range(9,21)]
    for i,row in enumerate(old['panels_data']):
        a,z=map(F,row['degrees']);assert a==F('1.50')+F(i,5000) and z==a+F(1,5000)
        c,num=newrig(a,z,w,w,endcap=F('.1'),min_cell=F(1,10**7))
        tot+=c
        assert (1-ll)*num[3]>3*ll and (1-ll)*num[4]>ll*(arb(3).sqrt()-arb(3)/2)
        tri.append(min(((1-ll)*num[3]-3*ll).lower(),((1-ll)*num[4]-ll*(arb(3).sqrt()-arb(3)/2)).lower()))
        s=b.degree(a).union(b.degree(z));mid=b.degree((a+z)/2);width=b.degree((z-a)/2)
        D=model.D_of(mid);Dj=model.D_of(arb_series([s,1],2))
        E3=g.line_pair(b.directions[4],g.unit(g.PI*7/6+s));center=(-E3[0],E3[1])
        ta=add(b.z,model.D_of(s));u0=sub(ta,center)
        # Q_A reflects along the 60 degree line.
        dy=(arb(3).sqrt()*u0[0]+u0[1])/2
        ex=g.HALF+(3*dy-4*dy**3)/16-dot(u0,ta)
        def harmonic_excess(par):
            e=g.line_pair(b.directions[4],g.unit(g.PI*7/6+par))
            p=add(b.z,model.D_of(par))
            u=sub(p,(-e[0],e[1]));dy=(arb(3).sqrt()*u[0]+u[1])/2
            return g.HALF+(3*dy-4*dy**3)/16-dot(u,p)
        ex=harmonic_excess(mid)+harmonic_excess(arb_series([s,1],2))[1]*(-width).union(width)
        assert ex>exact('.00273')
        harm.append(ex.lower())
        q0=sub(E3,mul(dot(E3,b.directions[1])-g.HALF,b.directions[1]))
        jQ=(-q0[0],q0[1]);Y=g.line_pair(b.directions[3],g.unit(g.PI*5/6+s))
        G=g.circle_line(Y,b.directions[5],1)
        P=sub(G,mul(dot(G,b.directions[2])-g.HALF,b.directions[2]))
        # The old polygon side P-F1 is a hexagon face. JQ is another
        # compulsory point of that face between P and F1.
        f1=g.line_pair(b.directions[2],b.directions[3])
        assert f1[0]<jQ[0]<P[0]
        arcang=arb.atan2(g.cross(ta,jQ),dot(ta,jQ))
        assert 0<arcang<exact('.28')
        m=[g.unit(g.PI/6+s+i*g.PI/3) for i in range(6)]
        mp=[g.unit(g.PI/6-s+i*g.PI/3) for i in range(6)]
        for v in ws:
            t=b.anchors(D,exact(v));tj=b.anchors(Dj,exact(v))
            pp=[tuple(p0+q[1]*(-width).union(width) for p0,q in zip(p,pj,strict=True)) for p,pj in zip(t,tj,strict=True)]
            gaps=[dot(m[0],pp[7]),dot(m[1],pp[6]),dot(mp[0],pp[7]),dot(mp[1],pp[6]),dot(mp[5],pp[0])]
            v1=(1-ll)*(min(x.lower() for x in gaps)-g.HALF)-ll/16
            assert v1>0
            cap.append(v1.lower())
        dm=exact(F(row['discrete_pose_lower']['exact'][0]))
        assert (1-ll)*dm>ll
        def tangent(left,right):
            par=b.degree(left).union(b.degree(right));d=model.D_of(par)
            raw=b.anchors(d,arb_series([exact(w),1],2))[6]
            q=tuple(v[0] for v in raw);v=tuple(t[1] for t in raw)
            out=mul(-1/g.norm(v),rot(v))
            val=dot(out,sub((arb(0),1/arb(3).sqrt()),q))
            if val>exact('.009'): return val.lower()
            assert right-left>F(1,10**8)
            center=(left+right)/2
            return min(tangent(left,center),tangent(center,right))
        sig.append(tangent(a,z))
    assert len(sig)==200
    delta0=ll*min(harm)
    deltaq=ll*min(sig)
    floor0=F(area['full_contour_area_minimum']['exact'][0]);rho=1-F('.268')*lam;floor=floor0-(1-rho**2)*F('.28')/6
    B=F('.844082046521218581802948589770767704354');assert floor>B
    return dict(complete=True,interval=old['interval'],panels=len(sig),lambda_exact=str(lam),q435_positive_root_boxes=tot,sector_radians_upper='.28',radial_retention=str(rho),
        root_margin_after_triangle=pack(min(tri)),five_caps_after_harmonic=pack(min(cap)),
        delta_0=pack(delta0.lower(),38),delta_q=pack(deltaq.lower(),38),
        area_lower=pack((floor,floor),45),excess_over_B=pack((floor-B,floor-B),44),
        source={p.name:sha256(p.read_bytes()).hexdigest() for p in (joint,strong,Path(inspect.getfile(rig)),Path(__file__))})
if __name__=='__main__':
    p=ArgumentParser()
    for k in ('comp','joint','strong','out'):p.add_argument(k,type=Path)
    a=p.parse_args();r=run(a.comp.resolve(),a.joint,a.strong);a.out.write_text(json.dumps(r,indent=1)+'\n')
    print(r['delta_0']['decimal'],r['delta_q']['decimal'],r['area_lower']['decimal'])
