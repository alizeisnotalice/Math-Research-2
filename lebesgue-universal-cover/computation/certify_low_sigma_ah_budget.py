"""Three Gibbs angles: forced straight-and-paired points cap the area saving of every anchored AH chord."""
import argparse
from fractions import Fraction as F
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys

import flint
from flint import arb_series,ctx
import certify_moving_sigma_fold_splice as model
from certify_low_sigma_fold_splice import ACTIVE,CYCLE
from certify_three_sigma_area_budget import exact,pack,subtract,digest,geometries,fixed_ah,ah_clearance,setup

b,g=model.b,model.g
add,sub,mul,dot,cross,rot,neg=b.add,b.sub,b.mul,b.dot,b.cross,b.rot,b.neg
CASES={
    '1.53': {'force':'.43042','fold':('.4303951507886','.4303951507888')},
    '1.52': {'force':'.43130','fold':('.4312779530696','.4312779530699')},
    '1.50': {'force':'.43310','fold':('.4330739295591','.4330739295595')},
}
CHAMBER=(F('52.65'),F('55.15'))
TURN=F('.019')
B=F('.844082046521218581802948589770767704354')
OTHER=tuple((k,s) for k in range(6) for s in (-1,1) if (k,s) not in ((0,1),(2,-1)))

def load(archive,name):
    return json.loads((archive/name).read_bytes())

def norm(p):
    return dot(p,p).sqrt()

def pull(u,k,s):
    c,v=g.unit(-k*g.PI/3)
    return (c*u[0]-v*u[1],s*(v*u[0]+c*u[1]))

def verify_fold(sigma,D,conf):
    left,right=map(F,conf['fold'])
    wforce=F(conf['force'])
    values=[model.delta(D,exact(w)) for w in (left,right,wforce)]
    assert values[0]>0 and values[1]<0 and values[2]<0
    derivative=model.delta(D,arb_series([b.ea(left,wforce),1],2))[1]
    assert derivative<0
    wbox=b.ea(left,wforce)
    points=b.anchors(D,wbox)
    angle=b.degree(CHAMBER[0]).union(b.degree(CHAMBER[1]))
    cone=[]
    for a,(ip,opposite) in enumerate(ACTIVE):
        i=CYCLE.index(ip)
        p=points[ip]
        dl,dr=[sub(p,points[CYCLE[(i+s)%len(CYCLE)]]) for s in (-1,1)]
        if cross(dl,dr)<0:
            dl,dr=dr,dl
        assert cross(dl,dr)>0
        u=g.unit(angle+2*a*g.PI/3)
        if opposite:
            u=neg(u)
        cone.extend([cross(dl,u),cross(u,dr)])
    assert all(value>0 for value in cone)
    v1,v2=model.harmonic(D,wbox)
    phi=g.arb.atan2(-v2,-v1)
    assert b.degree(CHAMBER[0])<phi<b.degree(CHAMBER[1])
    end_values=[v1*b.degree(edge).cos()+v2*b.degree(edge).sin()+g.HALF for edge in CHAMBER]
    assert all(value>0 for value in end_values)
    source=inspect.getsource(b.rig)
    token='scan(endcap,F(60)-endcap)'
    assert source.count(token)==1
    a,c=map(str,CHAMBER)
    source=source.replace(token,f'scan(endcap,F("{a}"));scan(F("{c}"),F(60)-endcap)')
    ns=b.rig.__globals__.copy()
    exec(compile(source,'rig_outside_low_sigma_fold','exec'),ns)
    parts,stats=ns['rig'](sigma,sigma,left,wforce,endcap=F('.1'),min_cell=F(1,1000000))
    assert min(stats)>0
    return dict(discriminant=tuple(map(pack,values)),negative_w_derivative=pack(-derivative),
                sector_cone_lower=pack(min(value.lower() for value in cone)),outside_angle_count=parts,
                chamber_endpoint_positive_minimum=pack(min(value.lower() for value in end_values)),
                end_slope_and_outside_lower=[pack(v) for v in stats])

def cap_wall(sigma,conf):
    s=b.degree(sigma)
    g.SIGMA=s
    n,m,*_=g.setup()
    D=model.D_of(s)
    L=norm(D)
    u=mul(-1/L,D)
    v=rot(u)
    p=F(conf['force'])
    w=exact(p)
    points=b.anchors(D,w)
    fold=verify_fold(sigma,D,conf)
    q=points[6]
    x=add(b.z,mul(w,D))
    q_series=b.anchors(D,arb_series([w,1],2))[6]
    qp=tuple(a[1] for a in q_series)
    nu=mul(-1/norm(qp),rot(qp))
    V=(g.arb(0),1/g.arb(3).sqrt())
    paired=dot(nu,sub(V,q))
    straight=-dot(v,sub(V,x))
    assert paired>0 and straight>0 and dot(u,qp)>0
    gaps=[]
    for k,s in OTHER:
        transformed=[pull(m[i],k,s) for i in (3,5)]
        gap=max((dot(z,a)-g.HALF).lower() for z in transformed for a in points)
        assert gap>0,(k,s,gap)
        gaps.append(gap)
    lam=dot(u,sub(q,b.z))/L
    eta=dot(v,sub(q,b.z))/L
    theta=(eta/(w+lam)).atan()
    assert lam>0 and eta>0 and 0<theta<exact(TURN)
    return dict(force_parameter=str(p),fold_bracket=conf['fold'],fold_witness=fold,
                paired_projection=pack(lam),paired_height=pack(eta),straight_projection=str(p),
                other_ten_violation_lower=pack(min(gaps)),paired_tangent_violation_per_triangle_weight=pack(paired),
                straight_tangent_violation_per_triangle_weight=pack(straight),
                maximal_profile_angle=pack(theta),strict_domination_angle=str(TURN))

def run(archive,sigma,out):
    ctx.prec=192
    conf=CASES[sigma]
    se=F(sigma)
    env={**os.environ,'LEBESGUE_SIGMA_DEGREES':f'{sigma},{sigma}'}
    shape=json.loads(subprocess.check_output([sys.executable,'-B',str(archive/'certify_splice_connectors.py')],env=env))
    assert shape['passed'] and shape['sigma_degrees_interval']==[str(se)]*2
    base=json.loads(subprocess.check_output([sys.executable,'-B',str(archive/'certify_gibbs_area_refined.py'),sigma,'10','48']))
    assert base['passed'] and F(base['sigma_degrees'])==se
    old=load(archive,'eh-two-line-cover-certificate.json')
    ehparms=load(archive,'eh-upper-chain-parameters.json')
    geo=geometries(sigma,old['parameters'],ehparms,g)
    eh=geo['new_eh']-geo['old_eh']
    assert len(geo['facets'])==11 and eh>0
    necessary=cap_wall(se,conf)
    p=F(conf['force'])
    gain,crossing,derivative=fixed_ah(archive,sigma,p,TURN,upper=F('.50'))
    ab=tuple(map(F,base['area_interval']))
    eh2=tuple(map(F,pack(eh)['exact']))
    bound=subtract(subtract(ab,gain),eh2)
    assert bound[0]>B
    loss=subtract(bound,(B,B))
    extreme=subtract(ab,gain)
    eh_budget=subtract(extreme,(B,B))
    cap_clearance=ah_clearance(
                        setup(sigma,g),
                        geo['points'][0],geo['points'][-1],g,p,TURN)
    assert cap_clearance>0
    sources=['certify_moving_sigma_fold_splice.py','certify_two_sigma_ideal_placement.py',
             'certify_low_sigma_fold_splice.py','certify_sigma_paired_rigidity.py','certify_three_sigma_area_budget.py',
             'certify_splice_connectors.py','certify_gibbs_area_refined.py','certify_ah_combined_area.py',
             'eh-two-line-cover-certificate.json','eh-upper-chain-parameters.json','audit_optimized_eh_chord.py',
             'certify_eh_dual_hinge.py']
    doc=dict(complete=True,sigma_degrees=sigma,threshold_B=str(B),necessary=necessary,
             dominating_cap_crossing=pack(crossing),crossing_slope_lower=str(derivative),
             AH_gain_upper=pack(gain),gibbs_area=pack(ab),eleven_EH_gain=pack(eh),
             dominating_AH_and_EH_clearance=pack(cap_clearance),
             area_lower_endpoint_enclosure_for_dominating_cut=pack(bound),
             area_gap_above_B=pack(loss),EH_gain_required_for_B_even_at_dominating_AH_cut=pack(eh_budget),
             shape_source=shape['script_sha256'],shape_receipt=shape,
             sources={name:digest(archive/name) for name in sources},
             script_sha256=digest(Path(__file__)),python=sys.version,flint=flint.__version__,precision=ctx.prec)
    out.write_text(json.dumps(doc,indent=2)+'\n')
    print(sigma,pack(bound)['decimal'],pack(loss)['decimal'],pack(eh_budget)['decimal'],
          necessary['maximal_profile_angle']['decimal'],necessary['fold_witness']['outside_angle_count'],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('archive',type=Path)
    ap.add_argument('sigma',choices=CASES)
    ap.add_argument('output',type=Path)
    a=ap.parse_args()
    run(a.archive,a.sigma,a.output)
