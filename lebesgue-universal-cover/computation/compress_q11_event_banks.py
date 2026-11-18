"""Coarsen frozen Q11 angular covers using exact edge events; keep all pose guards."""
from argparse import ArgumentParser
from fractions import Fraction as F
from pathlib import Path
import json
import sys

from flint import arb,ctx

from certify_q11_support_events import Events
from certify_mixed_support_budget import Host,Problem,ar,cloud,read,sha
from certify_high_six_host_neighborhood import Caps
from certify_mixed_leaf_centers import perimeter,verify_cap_radius,SCALE


def tupled(x):
    return tuple(map(tupled,x)) if isinstance(x,list) else x


def compress(old,evaluate):
    boundaries=[F(0)]
    keys=[]
    for row in old:
        i,k,*aux=row
        assert boundaries[-1]==F(i,2**k)
        boundaries.append(F(i+1,2**k))
        keys.append(tupled(aux))
    assert boundaries[-1]==1
    result=[]
    cursor=events=0
    minimum=arb(1)
    while cursor<len(keys):
        best=cursor
        chosen=None
        for key in dict.fromkeys(keys[cursor:cursor+12]):
            lo=best
            hi=len(keys)+1
            if lo>cursor:
                if boundaries[lo]-boundaries[cursor]>F(1,2):
                    continue
                if not evaluate(boundaries[cursor],boundaries[lo],key)[0]>0:
                    continue
            while hi-lo>1:
                m=(lo+hi)//2
                if boundaries[m]-boundaries[cursor]>F(1,2):
                    hi=m
                elif evaluate(boundaries[cursor],boundaries[m],key)[0]>0:
                    lo=m
                else:
                    hi=m
            if lo>best:
                best=lo
                chosen=key
        assert best>cursor and chosen is not None
        a,b=boundaries[cursor],boundaries[best]
        slack,n=evaluate(a,b,chosen)
        assert slack>0
        events+=n
        minimum=minimum.min(slack)
        result.append([str(a),str(b),*chosen])
        cursor=best
    return result,dict(original=len(old),compressed=len(result),events=events,minimum=str(minimum))


class MixedEvent:
    def __init__(self,host,ev,goal,radius):
        self.host,self.ev=host,ev
        self.correction=ar(host.area-goal-radius*(perimeter(host)/2+2*goal))
        self.lengths=[arb(s).sqrt()/host.den for s,_,_ in host.edges]
        self.base=host.offsets('mixed')
        self.r=ar(radius)
        self.models={}

    def __call__(self,l,r,key):
        # Unchanged frozen mixed weights; material indices are unused by the exact supports.
        packed,=tupled(key)
        if packed not in self.models:
            dirs,weighted=[],[]
            reciprocal=arb(0)
            used=set()
            for j,_,a,b in packed:
                w=F(a,b)
                assert 0<=j<len(self.base) and j not in used and 0<w<=F(1,2)
                used.add(j)
                dirs.append(self.host.normals[j])
                weighted.append(w)
                reciprocal+=ar(w)*(self.r*self.lengths[j]-self.base[j])
            self.models[packed]=dirs,weighted,reciprocal+self.correction
        dirs,weighted,offset=self.models[packed]
        lower,k=self.ev.minimum(l,r,dirs,weighted,self.host.den,False)
        return lower+offset,k


class CapEvent:
    def __init__(self,obj,ev):
        self.cap,self.ev=obj,ev

    def __call__(self,l,r,key):
        j,p,*_=key
        assert 0<=j<len(self.cap.triples)
        ids,weights,_=self.cap.triples[j]
        dirs=[self.cap.host.normals[i] for i in ids]
        lo,m=self.ev.minimum(l,r,dirs,weights,self.cap.host.den,False)
        return lo-self.cap.bound(p,j),m


def bank(comp,payload,coarsen=False):
    ctx.prec=160
    sys.path.insert(1,str(comp))
    Q=comp/'multi-host-points.json'
    raw=read(Q)
    assert payload['scale']==SCALE and sha(Q)==payload['Q11_sha256']
    pts,_=cloud(raw['points'],raw['point_denominator'])
    pi=arb.pi()
    goal=F(payload['target'])
    verify_cap_radius(goal,pi)
    summary={}
    for key,lot in payload['leaves'].items():
        assert len(lot['certs'])==len(lot['guards'])
        init_count=final_count=events=0
        old_minimum={}
        original_formula={}
        new_minimum={}
        for row,proof in zip(lot['guards'],lot['certs'],strict=True):
            coords,query=row[:-1],row[-1]
            assert len(coords)==11 and type(query) is int and query>0
            h=Host([F(z,SCALE) for z in coords],pi)
            assert h.hull_hash==proof['inner_hull_sha256']
            method=proof.get('strategy','mixed')
            if method=='cap':
                d=F(proof['radius'])
                assert F(query,SCALE)<=d
                cap=Caps(h,pts,pi,goal,d)
            else:
                assert method=='mixed'
                d=F(query,SCALE)
                if coarsen:
                    original=Problem(h,pts,'mixed',pi,F(0))
            for sign in (-1,1):
                ev=Events(raw['points'],raw['point_denominator'],sign,pi)
                fn=MixedEvent(h,ev,goal,d) if method=='mixed' else CapEvent(cap,ev)
                hand=str(sign)
                if coarsen:
                    source=proof['angles'][hand]
                    old_m=arb(1)
                    for i,k,*args in source:
                        l,r=F(i,2**k),F(i+1,2**k)
                        mv,_=fn(l,r,args)
                        if method=='cap':
                            earlier=cap.evaluate(l,r,sign,*args)
                        else:
                            earlier=fn.correction+original.evaluate(l,r,sign,args[0])
                            earlier+=sum((ar(d*F(a,b))*fn.lengths[j] for j,_,a,b in args[0]),arb(0))
                        assert earlier>0 and mv>0 and not mv<earlier
                        prev=original_formula.get(method,arb(1))
                        original_formula[method]=prev.min(earlier)
                        old_m=old_m.min(mv)
                    old_minimum[method]=old_minimum.get(method,arb(1)).min(old_m)
                    frozen,t=compress(source,fn)
                    proof['angles'][hand]=frozen
                    init_count+=t['original']
                else:
                    frozen=proof['angles'][hand]
                before=F(0)
                lower=arb(1)
                for l,r,*args in frozen:
                    l,r=F(l),F(r)
                    assert l==before<r<=1
                    v,m=fn(l,r,args)
                    assert v>0
                    lower=lower.min(v)
                    events+=m
                    before=r
                assert before==1
                new_minimum[method]=new_minimum.get(method,arb(1)).min(lower)
                final_count+=len(frozen)
        summary[key]=dict(old_count=init_count,new_count=final_count,events=events,
                          old_minima={k:str(x) for k,x in old_minimum.items()},
                          original_formula_minima={k:str(x) for k,x in original_formula.items()},
                          new_minima={k:str(x) for k,x in new_minimum.items()})
    return summary


def replay(comp,payload):
    assert payload['angle_verifier']=='events'
    bank(comp,payload)
    return payload


if __name__=='__main__':
    arg=ArgumentParser(description=__doc__)
    for name in ('comp','model','output'):
        arg.add_argument(name,type=Path)
    a=arg.parse_args()
    sys.path.insert(1,str(a.comp))
    data=read(a.model)
    d=bank(a.comp,data,True)
    data.update(angle_verifier='events',fixed_rows_source_sha256=sha(a.model),compression=d)
    a.output.write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(json.dumps(d,separators=(',',':')))
