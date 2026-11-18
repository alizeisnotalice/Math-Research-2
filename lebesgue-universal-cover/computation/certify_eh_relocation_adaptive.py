#!/usr/bin/env python3
"""Certify root-conditioned placement with a separate angular subdivision for each edge."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from certify_eh_relocation_root_condition import (
    RootConditionProgram, exact_fraction, phi_objective, root_clearance_objective,
)
from audit_eh_continuous_placement import polygon
from audit_optimized_eh_chord import PI, arb, minus, rot, scale, setup, trial, unit


def multipliers(label):
    base=("0","-1","-.5","1",".5","-1.5","-2","1.5","2")
    if label=="U1-Te":
        base=("0","-.2","-.175","-.25","-.3","-.35","-1","-.5","1",".5")
    return tuple(map(Fraction,base))


def certify(program,right,max_depth,replay=None):
    left=Fraction(30)
    objective=phi_objective(program,right,1)
    right_dual=program.optimize(objective) if replay is None else replay["right_end_dual"]
    sign_upper=program.validate(objective,right_dual)-arb(3)/2
    assert sign_upper<0, str(sign_upper)
    edges,_=polygon(program.delta,4)
    hex_names={"B1-Q","P-F1","F1-H","E3-D1","D1-G","C3-B1"}
    assert len([edge for edge in edges if edge[0] in hex_names])==6
    assert any(label=="K-C3" for label,_,_ in edges)
    _,_,c3,_,_,u,_,_,_=setup()
    co,si=unit(arb(program.delta))
    radial=minus(u,c3)
    u1_radial=minus(scale(co,radial),scale(si,rot(radial)))
    assert all(x<0 for x in u1_radial)
    results={}
    margins=[]
    for label,normal,cap in edges:
        # The exact third-quadrant C3 arc follows from the CB unit contact
        # after the chord/m5 and the original six hexagonal supports are certified.
        if label in hex_names or label.startswith("chord_H_"):
            continue
        cache={}
        used=set()
        cells=[]

        def endpoint(angle,mu):
            key=f"{angle}_{mu}"
            if key not in cache:
                objective,threshold=root_clearance_objective(program,angle,normal,cap,mu)
                dual=program.optimize(objective) if replay is None else replay["edges"][label]["duals"][key]
                cache[key]=(program.validate(objective,dual),threshold,dual)
            return key,cache[key]

        def margin(a,b,mu):
            ka,(u,threshold,_)=endpoint(a,mu)
            kb,(v,_,_)=endpoint(b,mu)
            cosine=(exact_fraction(b-a)*PI/360).cos()
            assert cosine>0 and max(u,v)>0
            bound=threshold-max(u,v)/cosine
            return bound,ka,kb

        def record(a,b,mu,w,ka,kb):
            assert w>0
            used.update((ka,kb))
            cells.append([str(a),str(b),str(mu)])
            margins.append((w,label))

        def split(a,b,depth):
            for mu in multipliers(label):
                w,ka,kb=margin(a,b,mu)
                if w>0:
                    record(a,b,mu,w,ka,kb)
                    return
            if depth==max_depth:
                raise ArithmeticError(f"cannot certify {label}: {a}--{b}")
            center=(a+b)/2
            split(a,center,depth+1)
            split(center,b,depth+1)

        start=len(margins)
        if replay is None:
            split(left,right,0)
        else:
            previous=left
            for aa,bb,mm in replay["edges"][label]["cells"]:
                a,b,mu=map(Fraction,(aa,bb,mm))
                assert a==previous and a<b<=right
                record(a,b,mu,*margin(a,b,mu))
                previous=b
            assert previous==right
        local_margin=min(w.lower() for w,_ in margins[start:])
        print(label,len(cells),str(local_margin),flush=True)
        results[label]={"cells":cells,"duals":{key:cache[key][2] for key in sorted(used)}}

    worst,label=min(margins,key=lambda row:row[0].lower())
    return {"summary":{"delta":program.delta,"angle_root_range_degrees":[str(left),str(right)],
                        "F_right_upper":str(sign_upper),"tested_nonhex_edges":len(results),
                        "C3_arc_retained_exactly":True,
                        "nonhex_slack_strict_lower":str(worst.lower()),"active_edge":label,
                        "angle_edge_cells":len(margins),"adaptive_cut_rows":len(program.adaptive_cuts),
                        "validated_duals":1+sum(len(r["duals"]) for r in results.values()),
                        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            "right_end_dual":right_dual,"adaptive_cuts":program.adaptive_cuts,"edges":results}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("delta")
    parser.add_argument("right_endpoint_degrees")
    parser.add_argument("output",type=Path)
    parser.add_argument("--max-depth",type=int,default=6)
    parser.add_argument("--replay",type=Path)
    args=parser.parse_args()
    geometry=trial(args.delta,tight=True)
    program=RootConditionProgram(args.delta)
    replay=None
    if args.replay:
        replay=json.loads(args.replay.read_text())
        assert replay["summary"]["delta"]==args.delta
        for i,g in replay["adaptive_cuts"]:
            program.cut(i,g,record=True)
    result=certify(program,Fraction(args.right_endpoint_degrees),args.max_depth,replay)
    result["summary"]["cut_gain_area_interval"]=geometry["gain_area"]
    args.output.write_text(json.dumps(result,separators=(",",":"))+"\n")
    print(json.dumps(result["summary"],indent=2),flush=True)
