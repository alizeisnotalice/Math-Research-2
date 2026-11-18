"""Freeze and replay mixed-area guards centered at exact pose-box midpoints."""
from argparse import ArgumentParser
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import json
import sys

from flint import arb,ctx

from certify_mixed_support_budget import Host,Problem,ar,cloud,read,sha

SCALE=2**45


def verify_cap_radius(target,pi):
    from certify_high_six_host_neighborhood import R
    r=ar(R)
    area=(pi-(1/(2*r)).acos())/4+(r*r-ar(F(1,4))).sqrt()/2
    assert ar(target)<area


def perimeter(host):
    poly=host.material
    count=0
    for a,b in zip(poly,poly[1:]+poly[:1],strict=True):
        sq=(a[0]-b[0])**2+(a[1]-b[1])**2
        x=isqrt(sq)
        count+=x+(x*x<sq)
    return F(count,host.den)


def verify_one(host,problem,rows,target):
    coeff=ar(perimeter(host)/2+2*target)
    lengths=[arb(sq).sqrt()/host.den for sq,_,_ in host.edges]
    least=arb(1)
    for branch in (-1,1):
        before=F(0)
        for i,k,terms in rows[str(branch)]:
            left,right=F(i,2**k),F(i+1,2**k)
            assert before==left<right<=1
            gain=problem.evaluate(left,right,branch,terms)
            assert gain>0
            mass=sum((ar(F(a,b))*lengths[j] for j,_,a,b in terms),arb(0))
            slope=coeff-mass
            assert slope>0
            least=least.min((ar(host.area-target)+gain)/slope)
            before=right
        assert before==1
    return least


def radius_rows(host,problem,target,radius):
    coeff=ar(perimeter(host)/2+2*target)
    lengths=[arb(sq).sqrt()/host.den for sq,_,_ in host.edges]
    def one(hand):
        def cell(i,k):
            left,right=F(i,2**k),F(i+1,2**k)
            row=problem.suggest(left,right,hand)
            if row:
                gain=problem.evaluate(left,right,hand,row)
                mass=sum((ar(F(a,b))*lengths[j] for j,_,a,b in row),arb(0))
                if gain>0 and ar(host.area-target)+gain>ar(radius)*(coeff-mass):
                    yield [i,k,row]
                    return
            assert k<15
            yield from cell(i*2,k+1)
            yield from cell(i*2+1,k+1)
        return [r for i in range(64) for r in cell(i,6)]
    return {str(hand):one(hand) for hand in (1,-1)}


def replay_bank(comp,payload):
    ctx.prec=160
    sys.path.insert(0,str(comp))
    Q11=comp/'multi-host-points.json'
    raw=read(Q11)
    assert payload['scale']==SCALE and sha(Q11)==payload['Q11_sha256']
    pts,_=cloud(raw['points'],raw['point_denominator'])
    p=arb.pi()
    goal=F(payload['target'])
    verify_cap_radius(goal,p)
    for key,lot in payload['leaves'].items():
        assert len(lot['certs'])==len(lot['guards'])
        for row,proof in zip(lot['guards'],lot['certs'],strict=True):
            coords,query=row[:-1],row[-1]
            assert len(coords)==11 and type(query) is int and query>0
            body=Host([F(z,SCALE) for z in coords],p)
            assert body.hull_hash==proof['inner_hull_sha256']
            if proof.get('strategy','mixed')=='mixed':
                low=verify_one(body,Problem(body,pts,'mixed',p,F(0)),proof['angles'],goal)
                assert ar(F(query,SCALE))<low
            else:
                from certify_high_six_host_neighborhood import Caps
                assert proof['strategy']=='cap'
                radius=F(proof['radius'])
                assert F(query,SCALE)<=radius
                obj=Caps(body,pts,p,goal,radius)
                for h in (1,-1):
                    obj.check(h,proof['angles'][str(h)])
    return payload


def make_bank(comp,proposal,target,wanted):
    ctx.prec=160
    sys.path.insert(0,str(comp))
    Q11=comp/'multi-host-points.json'
    raw=read(Q11)
    pts,_=cloud(raw['points'],raw['point_denominator'])
    pi=arb.pi()
    verify_cap_radius(target,pi)
    guards={}
    for key,items in proposal.items():
        rows=[]
        certs=[]
        for candidate in items:
            if len(rows)>=wanted:
                break
            values=candidate['center']
            body=Host([F(z,SCALE) for z in values],pi)
            strategy=candidate.get('strategy','mixed')
            proof=dict(inner_hull_sha256=body.hull_hash,sample=candidate.get('sample'),strategy=strategy)
            if strategy=='mixed':
                problem=Problem(body,pts,'mixed',pi,F('.000003'))
                try:
                    wanted_radius=F(candidate.get('proposed_radius','0'))
                    if wanted_radius:
                        side=radius_rows(body,problem,target,wanted_radius)
                    else:
                        side={str(h):problem.make(h) for h in (1,-1)}
                except AssertionError:
                    print('candidate angular certificate incomplete',key,candidate.get('sample'),file=sys.stderr,flush=True)
                    continue
                r=verify_one(body,problem,side,target)
                num=int(float(r.mid())*SCALE)-2
                threshold=F(candidate.get('required','0'))
                if num<=0 or ar(F(num,SCALE))<=ar(threshold):
                    continue
                assert ar(F(num,SCALE))<r
            else:
                from certify_high_six_host_neighborhood import Caps
                assert strategy=='cap'
                for q in (F(1),F(17,20),F(7,10)):
                    radius=F(candidate['radius'])*q
                    obj=Caps(body,pts,pi,target,radius)
                    try:
                        side={str(h):obj.make(h) for h in (1,-1)}
                        break
                    except AssertionError:
                        side=None
                if side is None:
                    continue
                proof['radius']=str(radius)
                num=int(radius*SCALE)
            rows.append([*values,num])
            proof['angles']=side
            certs.append(proof)
            print('accept',key,candidate.get('sample'),float(F(num,SCALE)),file=sys.stderr,flush=True)
        guards[key]=dict(guards=rows,certs=certs)
    return dict(target=str(target),scale=SCALE,Q11_sha256=sha(Q11),leaves=guards)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('comp',type=Path)
    p.add_argument('certificate',type=Path)
    p.add_argument('--generate',type=Path)
    p.add_argument('--target',type=F)
    p.add_argument('--guards',default=16,type=int)
    a=p.parse_args()
    if a.generate:
        out=make_bank(a.comp,read(a.generate),a.target,a.guards)
        a.certificate.write_text(json.dumps(out,separators=(',',':'))+'\n')
    bank=replay_bank(a.comp,read(a.certificate))
    print(json.dumps(dict(guard_counts={k:len(v['guards']) for k,v in bank['leaves'].items()},
                          certificate_sha256=sha(a.certificate)),separators=(',',':')))
