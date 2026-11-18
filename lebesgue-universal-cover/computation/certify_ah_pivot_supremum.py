"""Verify ideal-pose openness up to the compatibility fold, and a physical neighborhood for all shallow pivots."""
import argparse,json,sys,contextlib,io
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from flint import arb
from certify_ah_two_hit_budgets import geometry
from certify_ah_margin_envelope import replace_line


def run(archive,fold,out):
    f=json.loads(fold.read_bytes());old,src,_=geometry(archive,[(F('.427'),F('.0048'))],far=F('.183'),negative=F('.210'),eps=F('.0051'))
    early=json.loads((archive/'ah-fine-pivot-certificate.json').read_bytes())['bins'];assert [r['bin'] for r in early]==list(range(68));early=early[:65]
    rmin=min(arb(r['root_margin_lower']).lower() for r in early);pmin=min(arb(r['placement_margin_lower']).lower() for r in early);assert rmin>0 and pmin>0 and old['knots'][65]==F('.425')
    # All the actual physical chart gates hold on a fixed supercap for the limiting argument.
    for key,val in [('cuts=',"cuts=[(F('.429'),F('.0001'))]"),('wp, wf,',"wp, wf, wmax, eps = F('.429'),F('.040'),F('.43'),F('.001')"),('pbound = ',"pbound = F('.042')")]:src=replace_line(src,next(x for x in src.splitlines() if x.startswith(key)),val)
    tube={'__file__':old['__file__'],'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):exec(compile(src,tube['__file__'],'exec'),tube)
    ns=old;A=ns['A'];dot,perp,cross,rot,trig,sup=[ns[k] for k in ('dot','perp','cross','rotate','trig','supports')];normals=list(map(ns['av'],ns['d']['n']))
    rm=A(1);pm=A(1);dmax=-A(1);poses=0;worst=None
    for j in range(100):
        a=F('.425')+F(j,20000);b=a+F(1,20000);v=ns['centered_points'](a,b);e=[ns['sub'](v[i],v[k]) for i,k in ns['orders']]
        def coeffs(points):return [sum((-1)**k*dot(points[i],(perp if tang else lambda x:x)(normals[k])) for k,i in enumerate((2,3,0))) for tang in (False,True)]
        x,y=coeffs(v);dx,dy=coeffs(ns['point_jets'](ns['box'](a,b))[1]);der=2*(x*dx+y*dy);assert der<0;dmax=max(dmax,der.upper())
        phi=(y/x).atan()/ns['angle_unit'];assert A('5.5')<phi<A('6.4')
        for k,item in enumerate(((2,0),(3,1),(0,1))):
            u=rot(normals[0],trig(phi+A(60*k)));cc=[cross(edge,u) for edge in e];assert [z for m,z in enumerate(ns['sector']) if not cc[m]<0 and not cc[(m+1)%10]>0]==[item]
        left=(-ns['F_at'](v,e,F('5.0'))).lower();assert left>0;rm=min(rm,left)
        for i in range(28):
            lo=F(5)+F(i,20);hi=lo+F(1,20);ang=ns['box'](lo,hi);sn=trig(ang)
            h=[sup(v,e,sn,normal) for normal in normals[:2]]
            for label,u,c in ns['rows']:
                det=cross(normals[0],normals[1]);al,be=cross(u,normals[1])/det,cross(normals[0],u)/det
                m=(c-sup(v,e,sn,u)+al*(h[0]-A('.5'))+be*(h[1]-A('.5'))).lower();assert m>0,(j,i,label,m)
                if m<pm:pm,worst=m,[j,i,label]
            poses+=1
    receipt={'old_ideal_range':['0','.425'],'old_ideal_area_certificate_sha256':sha256((archive/'ah-fine-pivot-certificate.json').read_bytes()).hexdigest(),'old_ideal_root_and_pose_lower':[str(rmin),str(pmin)],'fold_certificate_sha256':sha256(fold.read_bytes()).hexdigest(),'new_section':['.425','.43'],'root_left_endpoint_degrees':'5','root_left_sign_lower':str(rm),'peak_range_degrees':['5.5','6.4'],'discriminant_derivative_upper':str(dmax),'full_angular_fit_window_degrees':['5','6.4'],'fit_cells':poses,'uniform_placement_lower':str(pm),'worst_placement_cell':worst,'physical_supercap':['.429','.0001'],'physical_profile_sha256':sha256(src.encode()).hexdigest(),'witness_source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'bits':128,'python':sys.version,'flint':ns['flint_version']}
    out.write_text(json.dumps(receipt,indent=2)+'\n');print(str(rm),str(pm),worst)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('fold',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();run(opt.archive.resolve(),opt.fold,opt.output)
