"""Exclude both handednesses and all rigid motions in the pinned-strip C153."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
import time

from flint import arb, ctx, __version__ as flint_version


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def geometry(comp):
    sys.path.insert(0,str(comp))
    from deformed_reuleaux_support import body,environment
    ns=environment(str(comp),'1.5497')
    ctx.prec=192
    return ns,body


def outer_rows(comp, config, ns, body, dirs):
    from certify_nonmonotone_ah_core import hull
    from certify_tapered_forced_support import rational
    core,pent,curves=[json.loads((comp/'finite-b-core'/n).read_bytes()) for n in ('core-outer.json','pentagons.json','curves.json')]
    assert pent['directions']==curves['N']
    q=core['vertex_grid']; assert q%pent['grid']==0 and pent['grid']==curves['scale']
    sc=q//pent['grid']
    pts=hull([tuple(p) for p in core['outer_integer_vertices']]+[(sc*x,sc*y) for x,y in pent['vertices_ccw_integer']+curves['combined_hull']])
    assert 2*max(F(pent['coordinate_rounding_error_strict']),F(curves['coordinate_error_asserted']))**2<F(3,10**11)**2
    delta=2*(arb.pi()/(2*pent['directions'])).sin()**2+rational(F(3,10**11))
    A=ns['A'];trig=ns['trig']; dot=ns['dot']; rt3=arb(3).sqrt()
    pins=(trig(A(0)),trig(A(120)))
    def calculate_placed_body(w,e,k,s):
        vs,_=body(ns,F(w),F(e))
        es=[ns['sub'](vs[i],vs[j]) for i,j in ns['orders']]
        def raw(u):
            return ns['supports'](vs,es,trig(A(-F(k))),(u[0],s*u[1]))
        a,b=(arb('.5')-raw(n) for n in pins)
        off=(a,(2*b+a)/rt3)
        return lambda n: raw(n)+dot(n,off)
    rowbody=[calculate_placed_body(*spec) for spec in config['bodies']]
    y,z=config['family_interval']['first'],config['family_interval']['last']
    # The pinned translations and the support function together have derivative Lipschitz constant 12.
    assert len(rowbody)==3 and (y,z)==(1,2) and config['translation_rule']=='pin_hexagon_faces_0_and_2'
    l,r=config['bodies'][y],config['bodies'][z]
    assert l[0]==r[0] and l[2:]==r[2:] and F(l[1])<F(r[1])
    pad=A(12)*A(F(r[1])-F(l[1]))**2/8
    out=[]
    for d in dirs:
        n=trig(A(F(d)))
        old=max((n[0]*x+n[1]*y).upper() for x,y in pts)/q+delta
        bds=[old.upper()]+[(fn(n)+(pad if i in (1,2) else 0)).upper() for i,fn in enumerate(rowbody)]
        top=max(bds)
        out.append((d,n,top))
    return out


def eliminate(ns, vs, es, rows, eta, thickness):
    A=ns['A'];trig=ns['trig'];box=ns['box'];supports=ns['supports']
    rt3=arb(3).sqrt(); half=arb('.5')
    pins=[trig(A(i)) for i in (0,60,120)]
    normrows=[(d,n,h,n[0]+n[1]/rt3,2*n[1]/rt3) for d,n,h in rows]
    def value(mid,b,s,derivative=False):
        return lambda n: supports(vs,es,b,(n[0],s*n[1]),derivative)
    def validate_angle(l,r,s):
        mid=(l+r)/2
        angle=trig(A(-mid)); rng=trig(box(-r,-l)); radius=A((r-l)/2)*ns['angle_unit']
        v=value(mid,angle,s); d=value(mid,rng,s,True)
        h0,h1,h2=map(v,pins);d0,d1,d2=map(d,pins)
        q=half+h1-h0-h2; dq=d0+d2-d1
        gap=radius*abs(dq)
        if q+gap<0: return 'strip-low',str((-(q+gap)).lower())
        if q-gap>A(eta):return 'strip-high',str((q-gap-A(eta)).lower())
        for label,n,z,a,b in normrows:
            pl=v(n)+a*(half-h0)+b*(half-h2)-z-A(thickness)
            slope=-d(n)+a*d0+b*d2
            residual=pl-radius*abs(slope)
            if residual>0:return label,str(residual.lower())
        return None,None
    data={};start=time.monotonic()
    for s in (1,-1):
        stack=[(F(i),F(i+1),0) for i in range(360)]
        count={};mins={};last=0
        while stack:
            l,r,depth=stack.pop()
            label,bound=validate_angle(l,r,s)
            if label is None:
                assert depth<26,(s,str(l),str(r),[(label,str(val))for label,n,val in rows])
                m=(l+r)/2
                stack.extend(((l,m,depth+1),(m,r,depth+1)))
                continue
            count[label]=count.get(label,0)+1
            val=arb(bound).lower();mins[label]=min(mins.get(label,val),val)
            last=max(last,depth)
        data[str(s)]=dict(leaves=sum(count.values()),max_depth=last,rows=count,residual_lower={k:str(v)for k,v in mins.items()})
        print('proved',s,data[str(s)]['leaves'],'leaves','depth',last, 'seconds',round(time.monotonic()-start,2),flush=True)
    return data


def main(a):
    comp=a.comp.resolve();ns,body=geometry(comp)
    config=json.loads(a.host.read_bytes()); trial=json.loads(a.input.read_bytes())
    rows=outer_rows(comp,config,ns,body,trial['directions'])
    values=[]
    for key in ('w','e'):
        middle=F(trial[key]);radius=F(trial.get(key+'_radius',0))
        assert radius>=0
        values.append(ns['box'](middle-radius,middle+radius))
    vs,_=body(ns,*values);es=[ns['sub'](vs[i],vs[j]) for i,j in ns['orders']]
    cert=eliminate(ns,vs,es,rows,a.eta,a.padding)
    result=dict(test={k:trial[k] for k in ('w','e','w_radius','e_radius') if k in trial},eta=str(a.eta),extra_support_padding=str(a.padding),
                rows={d:str(z)for d,n,z in rows},all_pose=cert,precision=ctx.prec,
                python=platform.python_version(),flint=flint_version,
                input_sha256=digest(a.input),host_sha256=digest(a.host),verifier_sha256=digest(Path(__file__)),
                base_sha256={n:digest(comp/'finite-b-core'/n) for n in ('core-outer.json','pentagons.json','curves.json')})
    a.out.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    if not __debug__:
        raise RuntimeError("Assertions must be enabled.")
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('comp',type=Path);a.add_argument('host',type=Path);a.add_argument('input',type=Path);a.add_argument('out',type=Path)
    a.add_argument('--eta',type=F,default=F('.000000301'));a.add_argument('--padding',type=F,default=F('.000003'))
    main(a.parse_args())
