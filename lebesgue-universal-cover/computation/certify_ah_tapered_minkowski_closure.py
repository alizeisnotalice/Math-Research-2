"""Prove three-pose coverage for the closed Minkowski convex hull of the R203 generators."""
from argparse import ArgumentParser
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import product
from pathlib import Path
import inspect
import json

from flint import arb, ctx
import certify_ah_tapered_release as tapered

LO=F(".4285")
HI=F(".95")
SCRIPT_SHA="ca131c3c44c07754f9abf7be22938e4d82d8c9942fc7b4c616f57fcd18b68733"
AUDIT_SHA="5658b6d51d60fbf36d15eaeba7d4b677716166a60aaccb9de9a2cda3d3e8f5d3"
POSES=[(240,-1),(0,1),(180,-1)]
SIMPLE_MARGIN=F("1e-9")
JOINT_MARGIN=F("8e-7")


def digest(raw):
    return sha256(raw).hexdigest()


def geometry(comp):
    source=comp/"certify_ah_tapered_release.py"
    assert digest(source.read_bytes())==SCRIPT_SHA
    assert digest(Path(tapered.__file__).read_bytes())==SCRIPT_SHA
    report=(comp/"ah-tapered-release.json").read_bytes()
    assert digest(report)==AUDIT_SHA
    old=json.loads(report)
    assert old["full_Minkowski_simplex"] and old["positions"]==65
    assert [tuple(p) for p in old["poses"]]==POSES
    assert tuple(map(F,old["w_interval"]))==(F(".427"),F(1))
    assert F(old["p"])==F(".95") and F(old["delta"])==F(".011")
    code=inspect.getsource(tapered.run).split("    def cell_slacks")[0]+"    return locals()\n"
    env=dict(tapered.run.__globals__)
    exec(compile(code,str(source),"exec"),env)
    g=env["run"](comp)
    assert g["poses"]==POSES and g["label"]=="taper_19"
    return g,digest(code.encode())


class Bounds:
    def __init__(self,prepared):
        self.ns=prepared["ns"]
        self.A=prepared["A"]
        self.g=prepared["g"]
        self.labels=[row[0] for row in prepared["tabs"][0]]
        assert len(self.labels)==len(set(self.labels))==189
        self.rows={}
        for sign,key in ((-1,"tabs"),(1,"tabs_positive")):
            self.rows[sign]={}
            for i,table in enumerate(prepared[key]):
                assert self.labels==[row[0] for row in table]
                self.rows[sign][i]={n:(u,h,tri) for n,u,h,tri in table}
        self.parts=((LO+(HI-LO)*F(j,32),LO+(HI-LO)*F(j+1,32)) for j in range(32))
        self.parts=tuple(self.parts)

    @lru_cache(maxsize=8192)
    def box(self,left,right):
        ns,g,A=self.ns,self.g,self.A
        m=(left+right)/2
        v=ns["points"](A(m))[0]
        vbox=ns["centered_points"](left,right)
        d=ns["point_jets"](A(left).union(A(right)))[1]
        edges=[g.minus(vbox[k],vbox[l]) for k,l in ns["orders"]]
        return v,d,edges

    @lru_cache(maxsize=200000)
    def terms(self,pose,label,left,right):
        ns,g=self.ns,self.g
        v,d,e=self.box(left,right)
        u,h,_=self.rows[-1][pose][label]
        tests=[g.cross(q,u) for q in e]
        selected=tuple(sorted({ns["sector"][j] for j in range(len(tests)) if
                                  not tests[j]<0 and not tests[(j+1)%len(tests)]>0}))
        assert selected
        return selected,tuple((g.dot(v[k],u)+c-h,g.dot(d[k],u)) for k,c in selected)

    def weighted_upper(self,terms,left,right):
        a=self.A
        choices=[]
        for i,label,c in terms:
            weight=a(c)
            _,possibly=self.terms(i,label,left,right)
            choices.append([(weight*x,weight*y) for x,y in possibly])
        rad=a((right-left)/2)
        return max((sum(x for x,_ in combo)+rad*abs(sum(y for _,y in combo))).upper()
                   for combo in product(*choices))

    def certify(self,terms,sign,threshold):
        a=self.A
        one=sum(a(c)*(self.rows[sign][i][label][2]-self.rows[sign][i][label][1]) for i,label,c in terms)
        assert one<-a(threshold),(terms,"triangle",str(one))
        todo=[(*row,0) for row in self.parts]
        leaves=0
        highest=one.upper()
        maxdepth=0
        while todo:
            left,right,depth=todo.pop()
            u=self.weighted_upper(terms,left,right)
            if u<-a(threshold):
                maxdepth=max(maxdepth,depth)
                leaves+=1
                highest=max(highest,u)
                continue
            assert depth<18,(terms,left,right,str(u),depth)
            middle=(left+right)/2
            todo.extend(((left,middle,depth+1),(middle,right,depth+1)))
        return dict(upper=str(highest),leaves=leaves,maxdepth=maxdepth)

    def identity_old_segment(self,prepared):
        # The x(w) vertex stays on the exact old Ta-Z supporting line.
        # The verified normal fan pins this same support vertex in every cell.
        ns,g,A=self.ns,self.g,self.A
        u,h,tri=self.rows[-1][1]["Ta-Z"]
        dep=g.dot(u,ns["z"])-h
        assert dep.contains(0) and g.dot(u,ns["D"]).contains(0)
        assert all((x-y).contains(0) for x,y in zip(u,g.neg(g.rot(ns["direct"]))))
        assert tri<h
        todo=[(*r,0) for r in self.parts]
        leaves=0
        while todo:
            l,r,k=todo.pop()
            if self.terms(1,"Ta-Z",l,r)[0]==((4,0),):
                leaves+=1
            else:
                assert k<18,(l,r,self.terms(1,"Ta-Z",l,r)[0])
                m=(l+r)/2
                todo.extend(((l,m,k+1),(m,r,k+1)))
        return {"leaves":leaves,"sector":[4,0],"triangle_strict_upper":str((tri-h).upper())}

    def antipodal(self):
        u,h,_=self.rows[-1][0]["K-C3"]
        v,k,_=self.rows[-1][2]["I-E3"]
        assert all((x+y).contains(0) for x,y in zip(u,v))
        assert (h-self.A(".5")).contains(0) and (k-self.A(".5")).contains(0)
        return {"rows":[[0,"K-C3"],[2,"I-E3"]],"weights":["1/2","1/2"],
                "weighted_support_violation":"0 by constant width one and antipodal preimages"}


def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("computation",type=Path)
    ap.add_argument("witness",type=Path)
    args=ap.parse_args()
    prep,src=geometry(args.computation)
    ctx.prec=160
    s=Bounds(prep)
    raw=args.witness.read_bytes()
    weights=json.loads(raw)["weights"]
    assert set(weights)==set(s.labels)-{"K-C3"}
    positive={n:s.certify([(0,n,F(1))],1,SIMPLE_MARGIN) for n in s.labels}
    identity=s.identity_old_segment(prep)
    neg_I={n:s.certify([(1,n,F(1))],-1,SIMPLE_MARGIN) for n in s.labels if n not in {"deep","Ta-Z"}}
    neg_alt={n:s.certify([(2,n,F(1))],-1,SIMPLE_MARGIN) for n in s.labels if n!="I-E3"}
    triples={}
    for label in s.labels:
        if label=="K-C3":
            continue
        lam=list(map(F,weights[label]))
        assert len(lam)==3 and min(lam)>=0 and sum(lam)==1
        terms=[(i,n,c) for i,(n,c) in enumerate(zip((label,"deep","I-E3"),lam,strict=True)) if c]
        triples[label]=s.certify(terms,-1,JOINT_MARGIN)
    print(json.dumps({
        "complete":True,"domain":[str(LO),str(HI)],"poses":POSES,
        "closed_Minkowski_convex_hull_certified":True,
        "source_script_sha256":SCRIPT_SHA,"source_audit_sha256":AUDIT_SHA,
        "derived_geometry_sha256":src,"witness_sha256":digest(raw),
        "positive_triangle_QA":positive,"negative_triangle_identity":neg_I,
        "negative_triangle_identity_Ta_Z":identity,"negative_triangle_alternative":neg_alt,
        "strict_Farkas_triples":triples,"exact_antipodal_triple":s.antipodal(),
        "strict_row_discharge_margin":str(SIMPLE_MARGIN),"strict_Farkas_margin":str(JOINT_MARGIN),
        "universal_placement_proved":False,
    },indent=2))


if __name__=="__main__":
    main()
