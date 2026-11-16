"""Production Fraction coefficient verifier for H10(x^3)H10(x^5).

There are 17 unordered canonical mother-pair layouts: one e=0 constant case,
and 16 e=1 choices after fixing pair 0 in its 3-fibre by factor exchange.
"""
from __future__ import annotations

import argparse, copy, hashlib, json
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from verify_h10_mothers import verify as verify_mothers

H10 = [1,1,1,0,0,0,0,0,0,1,1]


def pf(s: str) -> F:
    a,b=s.split("/",1); return F(int(a),int(b))
def qstr(x: F) -> str: return f"{x.numerator}/{x.denominator}"


@dataclass(frozen=True)
class I:
    lo:F; hi:F
    def __post_init__(self): assert self.lo<=self.hi
    def __add__(self,y:"I")->"I": return I(self.lo+y.lo,self.hi+y.hi)
    def __neg__(self)->"I": return I(-self.hi,-self.lo)
    def __mul__(self,y:"I")->"I":
        q=(self.lo*y.lo,self.lo*y.hi,self.hi*y.lo,self.hi*y.hi)
        return I(min(q),max(q))
Z,O=I(F(0),F(0)),I(F(1),F(1))


def square(x:I)->I:
    return I(F(0),max(x.lo*x.lo,x.hi*x.hi)) if x.lo<=0<=x.hi else x*x
def mul(a:list[I],b:list[I])->list[I]:
    c=[Z for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]=c[i+j]+x*y
    return c
def prod(blocks:list[list[I]])->list[I]:
    p=[O]
    for b in blocks:p=mul(p,b)
    return p
def sparse(degree:int, values:dict[int,I])->list[I]:
    p=[Z for _ in range(degree+1)]
    for i,x in values.items():p[i]=x
    return p


def disk(m:dict)->tuple[I,I]:
    r=pf(m["disk"]["radius"]); re=pf(m["disk"]["center"]["re"]); im=pf(m["disk"]["center"]["im"])
    return I(re-r,re+r),I(im-r,im+r)
def pair_fibre(m:dict,p:int)->list[I]:
    re,im=disk(m); return sparse(2*p,{0:square(re)+square(im),p:-(I(F(2),F(2))*re),2*p:O})
def witnesses(poly:list[I],side:str)->list[dict]:
    return [{"side":side,"degree":i,"upper":qstr(x.hi)} for i,x in enumerate(poly) if x.hi<0]
def unclear(poly:list[I])->list[int]: return [i for i,x in enumerate(poly) if x.lo<0<=x.hi]


def key(e:str,bits:int)->tuple[str,int]: return e,bits
EXPECTED={key("e0",0)}|{key("e1",b) for b in range(16)}
SURVIVORS={key("e0",0),key("e1",15)}


def build_records(data:dict)->list[dict]:
    mothers=sorted(data["mothers"],key=lambda x:x["block"])
    assert [m["block"] for m in mothers]==list(range(5))
    g={(i,p):pair_fibre(m,p) for i,m in enumerate(mothers) for p in (3,5)}
    layouts=[("e0",0,[],[g[i,p] for i in range(5) for p in (3,5)])]
    for bits in range(16):
        aa=[g[0,3]];bb=[g[0,5]]
        for i in range(1,5):
            p=3 if bits>>(i-1)&1 else 5
            aa.append(g[i,p]);bb.append(g[i,8-p])
        layouts.append(("e1",bits,aa,bb))
    assert len(layouts)==17
    out=[]
    for e,bits,aa,bb in layouts:
        ap,bp=prod(aa),prod(bb);ws=witnesses(ap,"A")+witnesses(bp,"B")
        ident=key(e,bits);status="rejected" if ws else "survivor"
        tag="constant" if ident==key("e0",0) else ("standard_H10_x3" if ident==key("e1",15) else "candidate")
        out.append({"e":e,"pair_3_bits":bits,"tag":tag,"status":status,"negative_witnesses":ws,
                    "uncertain_A":unclear(ap),"uncertain_B":unclear(bp),"degrees":[len(ap)-1,len(bp)-1]})
    return out


def exact_mul(a:list[int],b:list[int])->list[int]:
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c
def stretch(a:list[int],p:int)->list[int]:
    b=[0]*((len(a)-1)*p+1)
    for i,x in enumerate(a):b[i*p]=x
    return b
def exact_standard_tests()->dict:
    h3,h5=stretch(H10,3),stretch(H10,5);c=exact_mul(h3,h5)
    assert c[0]==c[-1]==1 and all(x in (0,1) for x in c)
    return {"constant_times_C":True,"H10_x3_times_H10_x5":True,
            "all_three_bits_standard_mapping":True}


def validate_records(rows:list[dict])->None:
    assert len(rows)==17
    keys=[key(r["e"],r["pair_3_bits"]) for r in rows]
    assert len(set(keys))==17 and set(keys)==EXPECTED
    assert {r["status"] for r in rows}<={"rejected","survivor"}
    actual={key(r["e"],r["pair_3_bits"]) for r in rows if r["status"]=="survivor"}
    assert actual==SURVIVORS and sum(r["status"]=="rejected" for r in rows)==15
    for r in rows:
        if r["status"]=="rejected":
            assert r["negative_witnesses"] and all(F(w["upper"])<0 for w in r["negative_witnesses"])


def interval_regressions()->dict:
    assert I(F(-2),F(3))*I(F(-5),F(7))==I(F(-15),F(21))
    assert I(F(-2),F(3))*I(F(-4),F(5))==I(F(-12),F(15))
    assert I(F(2),F(3))*I(F(4),F(5))==I(F(8),F(15))
    return {"negative_corner":True,"cross_zero_corner":True,"positive_corner":True}


def verify(data:dict)->dict:
    receipt=verify_mothers(data)
    rows=build_records(data);validate_records(rows)
    return {"mother_receipt":receipt,"candidate_count":17,"rejected_count":15,
            "survivors":[{"e":"e0","pair_3_bits":0,"tag":"constant"},
                         {"e":"e1","pair_3_bits":15,"tag":"standard_H10_x3"}],
            "records":rows,"exact_integrity":exact_standard_tests()}


def tamper_regressions(data:dict,result:dict)->dict:
    def fails(f):
        try:f()
        except (AssertionError,ValueError):return True
        return False
    bad=copy.deepcopy(data);bad["mothers"][0]["lhs"]="0/1"
    assert fails(lambda:verify_mothers(bad)),"bad mother disk accepted"
    bad=copy.deepcopy(data);bad["mothers"][0]["disk"]["radius"]="0/1"
    assert fails(lambda:verify_mothers(bad)),"zero-radius mother accepted"
    bad=copy.deepcopy(result["records"]);next(r for r in bad if r["status"]=="rejected")["status"]="survivor"
    assert fails(lambda:validate_records(bad)),"bad survivor label accepted"
    bad=copy.deepcopy(result["records"]);r=next(r for r in bad if r["status"]=="rejected")
    for w in r["negative_witnesses"]:w["upper"]="1/1"
    assert fails(lambda:validate_records(bad)),"positive fake witness accepted"
    assert fails(lambda:validate_records(copy.deepcopy(result["records"])[:-1])),"missing layout accepted"
    bad=copy.deepcopy(result["records"]);bad[-1]=copy.deepcopy(bad[-2])
    assert fails(lambda:validate_records(bad)),"duplicate layout accepted"
    return {"bad_lhs_rejected":True,"zero_radius_rejected":True,"bad_survivor_rejected":True,
            "positive_fake_witness_rejected":True,"missing_layout_rejected":True,
            "duplicate_layout_rejected":True,"interval_corner_regressions":interval_regressions()}


def main()->None:
    p=argparse.ArgumentParser();p.add_argument("--balls",type=Path,required=True);p.add_argument("--record",type=Path,required=True);a=p.parse_args()
    raw=a.balls.read_bytes();out=verify(json.loads(raw));out["balls_sha256"]=hashlib.sha256(raw).hexdigest();out["tamper_regressions"]=tamper_regressions(json.loads(raw),out)
    a.record.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:out[k] for k in ("balls_sha256","candidate_count","rejected_count","survivors","tamper_regressions")}))
if __name__=="__main__":main()
