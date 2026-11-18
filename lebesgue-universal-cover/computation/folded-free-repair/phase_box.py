"""Certify every remaining orientation of the earlier marked body under the free-area budget."""
import argparse, json, sys
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from flint import arb, arb_series, ctx, __version__ as flint_version

def run(comp,single,negative):
    sys.path.insert(0,str(single));from certify_counterexample import geometry
    from certify_monotone_cost import support_pins
    ns,body=geometry(comp)
    from certify_ah_tapered_minkowski_closure import geometry as get_geometry
    from certify_tapered_forced_support import rational as A
    ref=json.loads(negative.read_bytes())
    assert F(ref['area_budget'])==F('.844082046521218581802948589770767704354')
    assert F(ref['forbidden_negative_magnitude'])<=F('3.11e-7')
    old,_=get_geometry(comp);ctx.prec=192
    # R215 and the true area comparison give these extremes for any D with |D|<=B.
    # Check the conservative rounded constants directly from the exact stored face lengths.
    actual=[arb(z) for z in ref['old_supporting_segment_lengths']]
    U=F('.844082046521218581802948589770767704354')
    oldbase=F(json.loads((comp/'single-repair'/'single-repair.json').read_bytes())['base_area'][0])
    hostlo=F('.8440819208770037')
    upper=2*A(U-oldbase)/min(v.lower() for v in actual)
    bound=2*A(U-hostlo)/min(v.lower() for v in actual)
    assert upper<A(F('9.06e-7')) and bound<A(F('6.06e-7'))
    qa=A(F('9.06e-7')); qb=-A(F(ref['forbidden_negative_magnitude']))
    # n0,n2,n3,n5 are inside the original halfplanes on all of C*.
    delta=arb(0,A(F('6.06e-7')))
    g=old['g'];norms=[old['tilt'][j] for j in (3,5)]
    ic=json.loads((single/'cost-input.json').read_bytes())
    w0,e0=F(ic['w']),F(ic['e'])
    rw,re=F('0.0000004'),F('0.0000001')
    W,E=ns['box'](w0-rw,w0+rw),ns['box'](e0-re,e0+re)
    vb,_=body(ns,W,E);eb=[ns['sub'](vb[i],vb[j]) for i,j in ns['orders']]
    def jets(var):
        w=arb_series([W,int(var=='w')]);e=arb_series([E,int(var=='e')])
        l,t=ns['chart'](ns['add'](ns['z'],ns['mul'](w,ns['D'])))
        p=ns['psi'](l);t-=e
        add,mul,R=ns['add'],ns['mul'],ns['reflection']
        pts=[add(ns['l0'],mul(l,ns['el'])),R(add(ns['l0'],mul(p,ns['el']))),
             R(ns['X'](p,-t)),add(ns['mid'],mul(t,ns['en'])),ns['X'](l,t)]
        return [(a[1],b[1]) for a,b in pts]
    dw,de=jets('w'),jets('e')
    vs,_=body(ns,w0,e0);es=[ns['sub'](vs[i],vs[j]) for i,j in ns['orders']]
    roots=[ns['trig'](ns['A'](alpha)) for alpha in (0,60,120)]
    regions={'1':[['0','.00085'],['354.01','354.28'],['359.9975','360']],'-1':[['239.9975','240.00085']]}
    counts={}
    for s in (1,-1):
        allowed=[tuple(map(F,p)) for p in regions[str(s)]]
        edges=sorted({F(k) for k in range(361)}|{a for p in allowed for a in p})
        stack=[(a,b,0) for a,b in zip(edges,edges[1:])];stats={};depth=0
        while stack:
            l,r,k=stack.pop();mid=(l+r)/2;key=None
            if any(i<=l and r<=j for i,j in allowed):key='permitted_interval'
            if key is None:
                om=ns['box'](l,r);halfwidth=arb(0,ns['A']((r-l)/2)*ns['angle_unit'])
                fn,t,q0=support_pins(ns,vs,es,s,ns['A'](mid));trig=ns['trig'](-om)
                def deriv(n):return -ns['supports'](vb,eb,trig,(n[0],s*n[1]),True)
                def param(n):
                    co,si=ns['trig'](ns['A'](mid))
                    vn=(co*n[0]+s*si*n[1],-si*n[0]+s*co*n[1])
                    votes=[ns['cross'](p,vn) for p in eb]
                    ids=[ns['sector'][j][0] for j,h in enumerate(votes) if not h<0 and not votes[(j+1)%10]>0]
                    assert ids
                    def family(jets):
                        vals=[ns['dot'](vn,jets[i]) for i in ids];a=vals[0]
                        for b in vals[1:]:a=a.union(b)
                        return a
                    return family(dw),family(de)
                pp=list(map(param,roots));radii=[arb(0,A(rw)),arb(0,A(re))]
                a,b,c=map(deriv,roots);q=q0+halfwidth*(b-a-c)+sum(rad*(pp[1][j]-pp[0][j]-pp[2][j]) for j,rad in enumerate(radii))
                if q>qa:key='positive_face_budget'
                elif q<qb:key='mixed_negative_budget'
                else:
                    td=(-a,(-a-2*c)/arb(3).sqrt())
                    for n in norms:
                        gap=fn(n)+g.dot(n,t)-arb('.5')+halfwidth*(deriv(n)+g.dot(n,td))
                        aw=n[0]+n[1]/arb(3).sqrt();bw=2*n[1]/arb(3).sqrt()
                        pd=param(n);gap+=sum(rad*(pd[j]-aw*pp[0][j]-bw*pp[2][j]) for j,rad in enumerate(radii))
                        gap+=delta*aw+delta*bw
                        if gap*arb('.2535')/2>A(U-hostlo):key='oblique_area';break
            if key is None:
                assert k<29,(s,str(l),str(r),str(q))
                stack.extend(((l,mid,k+1),(mid,r,k+1)));continue
            stats[key]=stats.get(key,0)+1;depth=max(depth,k)
        counts[str(s)]=dict(leaves=stats,max_depth=depth)
    return dict(parameter_radius=[str(rw),str(re)],all_angles=True,regions_degrees=regions,classes=counts,witness=ic,
                extra_body=ref['witness'],q_negative_lower=ref['forbidden_negative_magnitude'],
                q_upper=str(upper),xy_budget=str(bound),arb_precision=ctx.prec,flint=flint_version,
                dependencies={p.name:sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),negative,single/'cost-input.json')})

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('comp','single','negative','out'):p.add_argument(name,type=Path)
    a=p.parse_args();result=run(a.comp.resolve(),a.single.resolve(),a.negative)
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(result['classes'])
