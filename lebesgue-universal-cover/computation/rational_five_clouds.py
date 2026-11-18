"""Build explicit diameter-one rational test sets and transfer the full five-body bound."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import ceil
from pathlib import Path
from certify_nine import S, append_point, ar, arb, arc_body, ctx

HERE=Path(__file__).resolve().parent
TARGET=F(120653235457106537,144115188075855872)
PI=F(22,7)

def build(disk_count,arc_count,goal):
    clouds=[]
    for n in (0,3,5,7,9):
        pts=[]
        if n:
            arc_body(pts,ceil(F(arc_count,n)),arb.pi(),n,arb(0),arb(0),arb(0))
        else:
            for j in range(disk_count):
                s,c=(2*arb.pi()*j/disk_count).sin_cos()
                append_point(pts,c/2,s/2)
        pts=sorted(set(pts))
        diameter=max((x-a)**2+(y-b)**2 for (x,y),(a,b) in combinations(pts,2))
        assert diameter<=(S+10)**2
        clouds.append(dict(n=n,points=pts,maximum_squared_distance_numerator=diameter))
    e=max(PI**2/(8*arc_count**2),PI**2/(4*disk_count**2))+F(12,S)
    bound=(TARGET-PI*e**2)/(1+2*e/(F(1,2)-e))
    assert bound>F(goal)
    obj=dict(denominator=S+10,clouds=clouds)
    raw=(json.dumps(obj,separators=(',',':'))+'\n').encode()
    return raw,dict(disk=disk_count,arc_resolution=arc_count,counts=[len(c["points"]) for c in clouds],
                    error=str(e),area_lower=str(bound),strict_decimal=goal,sha256=sha256(raw).hexdigest())

if __name__=="__main__":
    ctx.prec=160
    assert arb.pi()<ar(PI)
    parent=json.loads((HERE/"tomar-grouped-08372-public-audit.json").read_bytes())
    assert parent["complete"] and parent["target"]=="0.8372" and TARGET>F(parent["target"])
    summary=[]
    for a,b,goal in [(1024,648,".83719015"),(288,204,".83710030")]:
        raw,row=build(a,b,goal)
        filename=f"five-rational-clouds-{sum(row['counts'])}.json"
        if (HERE/filename).exists():
            assert (HERE/filename).read_bytes()==raw
        else:
            (HERE/filename).write_bytes(raw)
        summary.append(row)
    print(json.dumps(summary,indent=2))
