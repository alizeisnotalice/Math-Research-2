"""Integrate the compulsory E-H-T chain in the original Gibbs host."""
from argparse import ArgumentParser
from fractions import Fraction as F
import json
from hashlib import sha256
from pathlib import Path
import sys
from flint import arb, arb_series, ctx

def run(comp, joint, tips):
    sys.path.insert(0,str(comp))
    from certify_low_sigma_ah_budget import b,g,sub,mul,dot,cross
    from certify_three_sigma_area_budget import exact,pack
    ctx.prec=192
    ctx.cap=8
    def gain(sig):
        n=b.directions
        m=[g.unit(g.PI/6+sig+i*g.PI/3) for i in range(6)]
        c=g.line_pair(n[0],m[5])
        h=(-c[0],c[1])
        b3=g.line_pair(n[1],m[0])
        r=g.line_pair(m[2],n[1],g.HALF,0)
        s=g.circle_line(r,n[0],-1)
        e=sub(s,n[0])
        t=g.circle_pair(s,b3,1)
        ve,vt=sub(e,s),sub(t,s)
        cr=cross(ve,vt)
        dt=dot(ve,vt)
        at=lambda z:z[0] if isinstance(z,arb_series) else z
        assert at(dt)>0 and at(cr)>0
        arc=cross(s,sub(t,e))+(cr/dt).atan()
        return (cross(e,h)+cross(h,t)-arc)/2,cross(e,h),cross(h,t)
    prev=json.loads(joint.read_bytes())
    certified=json.loads(tips.read_bytes())
    assert prev['interval']==['3/2','77/50'] and prev['complete'] and certified['complete']
    stages=prev['panels_data']
    assert len(stages)==200
    derivatives=[]
    radial=[]
    forced_derivatives=[]
    for i,r in enumerate(stages):
        l,z=map(F,r['degrees'])
        assert l==F('1.50')+F(i,5000) and z==l+F(1,5000)
        a,mid,end=[b.degree(w) for w in (l,(l+z)/2,z)]
        slab=a.union(end)
        v,eh,ht=gain(arb_series([mid,1],5))
        rem,reh,rht=gain(arb_series([slab,1],5))
        delta=(a-mid).union(end-mid)
        for germ,whole in [(eh,reh),(ht,rht)]:
            crossing=germ[0]+germ[1]*delta+germ[2]*delta**2+whole[3]*delta**3
            assert crossing>0
            radial.append(crossing.lower())
        d=v[1]+2*v[2]*delta+3*v[3]*delta**2+4*rem[4]*delta**3
        derivatives.append(d)
        lower,upper=map(F,r['area_derivatives_per_radian'][0]['exact'])
        full=d+exact(upper)
        assert full<0
        forced_derivatives.append(full.upper())
    val=gain(b.degree(F('1.54')))[0]
    priorlo,priorhi=map(F,prev['universal_compulsory_area_minimum']['exact'])
    new=exact(priorlo).union(exact(priorhi))+val
    threshold=exact('.844082046521218581802948589770767704354')
    assert new-threshold>exact('.000002')
    return dict(complete=True,slabs=200,interval=prev['interval'],gain_at_minimum=pack(val,48),
       full_contour_area_minimum=pack(new,48),minimum_minus_B=pack(new-threshold,48),
       radial_cross_lower=pack(min(radial),35),maximum_area_derivative_upper=pack(max(forced_derivatives),35),
       dependencies={path.name:sha256(path.read_bytes()).hexdigest() for path in (joint,tips)},
       verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=='__main__':
    p=ArgumentParser()
    for k in ('comp','joint','tips','output'):p.add_argument(k,type=Path)
    a=p.parse_args()
    r=run(a.comp.resolve(),a.joint,a.tips)
    a.output.write_text(json.dumps(r,indent=1)+'\n')
    print(r['full_contour_area_minimum']['decimal'],r['minimum_minus_B']['decimal'])
