"""Extend the C1--C4 tube to C2 after exact simultaneous dihedral alignment."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
import json
from pathlib import Path

from flint import arb, ctx

import certify_six_anchor_path as first
from certify_mixed_support_budget import PI_UP, ar, cloud

CENTERS=32
AREA_PIECES=256
CURVED_SPEED=F(".00127")
P30_SPEED=F(".00192")
RATIONAL_ERROR=F(1,10**12)
TAIL_ERROR=2*RATIONAL_ERROR
AREA_SHA="39e878905779fa8fc54bc18e4739475b55c1e40eb0a18a6f2791d25665f941ef"
Q11_SHA="dd524963c0a92dadeee96ed5fc511bb85db1c52ec6ee672258038591b97c61e4"

def configure():
    first.CENTERS=CENTERS
    first.AREA_PIECES=AREA_PIECES
    first.FIVE_SPEED=CURVED_SPEED
    first.POINT_SPEED=P30_SPEED

def transformed(pose,pi):
    v=list(map(ar,pose))
    result=v.copy()
    rt=arb(3).sqrt()/2
    for n,j,k in ((3,-1,0),(5,2,3),(7,5,3),(9,8,4),(1,11,0)):
        if j>=0:
            result[j]=4*pi/3-v[j]-2*k*pi/n
        x,y=v[j+1:j+3]
        result[j+1:j+3]=[-x/2-rt*y,-rt*x+y/2]
    return result

def rational_endpoint(old,pi,stored):
    x=transformed(old,pi)
    if stored is None:
        numbers=[F(round(float(a.mid())*10**12),10**12) for a in x]
    else:
        numbers=list(map(F,stored))
    assert all(abs(a-ar(b))<ar(RATIONAL_ERROR) for a,b in zip(x,numbers,strict=True))
    return numbers

def tail(comp):
    # The exact aligned endpoint is isometric to C4, so both published
    # estimates transfer. The rational endpoint correction costs 2e-12.
    refs=[]
    for name,pinned in (("structural-six-escape-audit.json",AREA_SHA),
                        ("q11-structural-anchor-audit.json",Q11_SHA)):
        data=(comp/name).read_bytes()
        assert first.sha(data)==pinned
        refs.append(json.loads(data))
    xyz,cap=refs
    raw=(comp/"structural-six-escape-input.json").read_bytes()
    assert xyz["input_sha256"]==cap["config"]["inputs"]["anchor"]==first.sha(raw)
    assert arb(2).sqrt()+ar(F(".533"))<2
    lo,hi=map(F,xyz["area_bracket"])
    assert hi<first.TARGET
    d=F(cap["config"]["augmented_radius"])
    assert d==F(".00035") and F(cap["config"]["augmented_target"])==first.TARGET
    low=(lo-PI_UP*TAIL_ERROR**2)/(1+4*TAIL_ERROR)
    high=hi*(1+4*TAIL_ERROR)+PI_UP*TAIL_ERROR**2
    assert first.BASE<low<high<first.TARGET
    return {"maximum_Hausdorff_excursion":str(TAIL_ERROR),
            "radius":str(d-TAIL_ERROR),"area_bracket":[str(low),str(high)]}

def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("comp",type=Path)
    cli.add_argument("witness",type=Path)
    cli.add_argument("--generate",action="store_true")
    cli.add_argument("--workers",type=int,default=1)
    a=cli.parse_args()
    assert 1<=a.workers<=16
    configure()
    ctx.prec=160
    pi=arb.pi()
    assert pi<ar(PI_UP)
    data=first.data_files(a.comp)
    new=data["structural-six-escape-input.json"]
    assert new["hand"]==1
    old=list(map(F,new["pose"]))
    start=list(map(F,data["free-p30-placement.json"]))
    p30=data["finite-point-leaf-input.json"]
    p=set(map(tuple,p30["points"]))
    # The global reflection changes the formal hand, but the exact point
    # set is unchanged under reflection in its own x axis.
    assert p=={(x,-y) for x,y in p}
    Q=data["multi-host-points.json"]
    q11,d2=cloud(Q["points"],Q["point_denominator"])
    assert len(q11)==11
    conf={"inputs":first.EXPECTED,"extra_audits":{"C4_area":AREA_SHA,"C4_Q11":Q11_SHA},
          "alignment":"reflect in x axis then rotate by 4pi/3; angles of B5,B7,B9 shifted by 3,3,4 periods",
          "formal_P30_hand":1,"caps_centers":CENTERS,"area_cells":AREA_PIECES,
          "five_speed":str(CURVED_SPEED),"P30_speed":str(P30_SPEED),
          "base_area":str(first.BASE),"target":str(first.TARGET)}
    check=None if a.generate else json.loads(a.witness.read_bytes())
    if check:
        assert check["config"]==conf
    end=rational_endpoint(old,pi,check["end"] if check else None)
    speed,p_speed=first.motion(start,end,p30,pi)
    r=ar(F(1,2))
    R=ar(F(".703"))
    assert (pi-(r/R).acos())*r**2+r*(R**2-r**2).sqrt()>ar(first.TARGET)
    old_rows=[None]*CENTERS if check is None else check["centers"]
    assert len(old_rows)==CENTERS
    cc=[(first.midpoint(i,CENTERS),start,end,p30,q11,w) for i,w in enumerate(old_rows)]
    ac=[(i,AREA_PIECES,start,end) for i in range(AREA_PIECES)]
    if a.workers==1:
        cdata=list(map(first.capsule_job,cc)); areas=list(map(first.area_job,ac))
    else:
        with ProcessPoolExecutor(max_workers=a.workers,initializer=configure) as pool:
            cdata=list(pool.map(first.capsule_job,cc))
            areas=list(pool.map(first.area_job,ac,chunksize=4))
    if a.generate:
        raw=json.dumps({"config":conf,"end":list(map(str,end)),"centers":[w for _,w in cdata]},
                       separators=(",",":"))+"\n"
        tmp=a.witness.with_suffix(".pending")
        tmp.write_text(raw)
        tmp.replace(a.witness)
    suffix=tail(a.comp)
    c=[v for v,_ in cdata]
    rad=min(F(suffix["radius"]),*(F(z["interval_radius"]) for z in c))
    lo=min(F(suffix["area_bracket"][0]),*(F(x["whole_area_lower"]) for x in areas))
    hi=max(F(suffix["area_bracket"][1]),*(F(x["whole_area_upper"]) for x in areas))
    print(json.dumps({"complete":True,"config":conf,"witness_sha256":first.sha(a.witness.read_bytes()),
            "rational_endpoint":list(map(str,end)),"Q11_squared_distance_numerator":d2,
            "body_speeds":speed,"P30_speed":p_speed,
            "uniform_area_bracket":[str(lo),str(hi)],"uniform_area_display":[float(lo),float(hi)],
            "minimum_tube_radius":str(rad),"minimum_tube_radius_display":float(rad),
            "minimum_segment_P30_margin":min(d["p30_interval_slack_display"] for d in c),
            "tail":suffix,"tube_cells":c,"area_cells":areas},indent=2))

if __name__=="__main__":
    main()
