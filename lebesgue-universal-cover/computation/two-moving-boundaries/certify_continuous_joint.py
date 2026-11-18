"""Certify all Gibbs angles throughout a closed interval for the joint AH/EH contour."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from argparse import ArgumentParser
import json,sys
from flint import arb,arb_series,ctx
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here))
from continuous_curves import geometry

def run(comp):
    sys.path.insert(0,str(comp))
    from certify_low_sigma_ah_budget import model,exact,b,g,add,sub,mul,dot,cross,rot,OTHER,pull
    from certify_sigma_mandatory_arcs import anchors,certify_cell
    from certify_sigma_paired_rigidity import rig
    from certify_three_sigma_area_budget import pack
    ctx.prec=192;ctx.cap=8
    lo,hi=F('1.50'),F('1.54')
    ws=[F('.435')]+[F(i,20) for i in range(9,21)]
    count=200;rows=[]
    for part in range(count):
        a=lo+(hi-lo)*part/count;z=lo+(hi-lo)*(part+1)/count
        l,m,r=[b.degree(v) for v in (a,(a+z)/2,z)]
        sig=l.union(r);step=(l-m).union(r-m)
        arc_c=[];arc_r=[]
        central=geometry(arb_series([m,1],5),arc_c);outer=geometry(arb_series([sig,1],5),arc_r)
        for center,full in zip(arc_c,arc_r,strict=True):
            degree=3
            cr=center[0]+sum(center[k]*step**k for k in range(1,degree))+full[degree]*step**degree
            assert cr>0
        vds=[]
        for v,rem in zip(central,outer,strict=True):
            dv=v[1]+2*v[2]*step+3*v[3]*step**2+4*rem[4]*step**3
            assert dv<0
            vds.append(dv)
        nums=[];margins=[]
        for w in ws:
            num,st=rig(a,z,w,w,endcap=F('.1'),min_cell=F(1,10**6))
            assert min(st)>0
            nums.append(num);margins.append(min(st))
        ncell,ti=certify_cell('TI',a,z)
        gg=anchors('TI')[0];E,S,T,B=(gg[k] for k in ('E','S','T','B'))
        va,vb=sub(E,S),sub(T,S)
        assert dot(va,vb)>0 and cross(va,vb)>0
        angle=arb.atan2(cross(va,vb),dot(va,vb))
        assert 0<angle<exact('.04')
        cc,ss=g.unit(arb(0).union(angle))
        radial=add(mul(cc,va),mul(ss,rot(va)))
        dB=2*dot(sub(S,B),rot(radial))
        assert dB>0 and ti['lens_distance_slack']>0 and cross(vb,sub(T,B))>0
        mid=add(S,mul(1/g.norm(add(va,vb)),add(va,vb)))
        chord=sub(T,E);out=mul(-1/g.norm(chord),rot(chord))
        squared=dot(chord,chord)/4
        sagitta=squared/(1+(1-squared).sqrt())
        signed=sagitta-exact('.0000001')*(g.norm(S)+1)
        assert signed>exact('.00012')
        D=model.D_of(sig);pp=b.anchors(D,exact(ws[0]))
        raw=b.anchors(D,arb_series([exact(ws[0]),1],2))[6]
        vel=tuple(v[1] for v in raw);nu=mul(-1/g.norm(vel),rot(vel))
        straight=add(b.z,mul(exact(ws[0]),D))
        vert=(arb(0),1/arb(3).sqrt())
        pose=[dot(nu,sub(vert,pp[6])),-dot(rot(mul(-1/g.norm(D),D)),sub(vert,straight))]
        _,norms,*_=g.setup()
        for t,hand in OTHER:
            outward=[pull(norms[i],t,hand) for i in (3,5)]
            pose.append(max((dot(v,p)-g.HALF).lower() for v in outward for p in pp))
        assert min(pose[1:])>0
        def tangent_lower(left,right):
            ds=model.D_of(b.degree(left).union(b.degree(right)))
            qj=b.anchors(ds,arb_series([exact(ws[0]),1],2))[6]
            v=tuple(z[1] for z in qj);q=tuple(z[0] for z in qj)
            value=-dot(rot(v),sub(vert,q))/g.norm(v)
            if value>0:return value.lower()
            assert right-left>F(1,10**8)
            mm=(left+right)/2
            return min(tangent_lower(left,mm),tangent_lower(mm,right))
        pose[0]=tangent_lower(a,z)
        assert min(pose)>0
        rows.append(dict(degrees=[str(a),str(z)],root_boxes=sum(nums),TI_boxes=ncell,
                         root_margin_lower=pack(min(margins).lower()),TI_margin_lower=pack(min(ti[k] for k in ti if 'negative' not in k).lower()),
                         discrete_pose_lower=pack(min(pose).lower()),ET_witness_distance_lower=pack(signed.lower()),
                         area_derivatives_per_radian=[pack(q) for q in vds]))
    exact_min=geometry(b.degree(hi))[0];failed_max=geometry(b.degree(lo))[1]
    threshold=exact('.844082046521218581802948589770767704354')
    assert exact_min>threshold>failed_max
    result=dict(complete=True,interval=[str(lo),str(hi)],panels=count,closed_root_boxes=sum(r['root_boxes'] for r in rows),
                closed_TI_boxes=sum(r['TI_boxes'] for r in rows),
                universal_compulsory_area_minimum=pack(exact_min,36),minimum_minus_B=pack(exact_min-threshold,36),
                rejected_chord_area_maximum=pack(failed_max,36),panels_data=rows,
                source_sha256={name:sha256((here/name).read_bytes()).hexdigest() for name in ('certify_continuous_joint.py','continuous_curves.py')})
    return result

if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('computation',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    result=run(a.computation.resolve())
    a.output.write_text(json.dumps(result,indent=1)+'\n')
    print(result['closed_root_boxes'],result['universal_compulsory_area_minimum']['decimal'],flush=True)
