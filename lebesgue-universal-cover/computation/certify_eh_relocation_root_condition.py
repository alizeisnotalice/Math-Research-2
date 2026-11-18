#!/usr/bin/env python3
"""Probe a second relocation with the compatibility root retained in the clearance objectives."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np

from certify_eh_general_relocation import (
    Program, exact_fraction, numeric, phi_objective, rotate_angle,
)
from audit_eh_continuous_placement import polygon
from audit_optimized_eh_chord import HALF, PI, arb, ctx, dot, minus, plus, rot, scale, setup, trial, unit


class RootConditionProgram(Program):
    def __init__(self,delta,avoid_ah=True):
        super().__init__(delta)
        previous=self.dimension
        self.dimension+=2
        eye=np.eye(self.dimension+1,dtype=object)

        def expand(point):
            return np.column_stack((point[:,:previous],np.zeros((2,2),dtype=object),point[:,-1]))

        self.witnesses=[expand(p) for p in self.witnesses]
        self.support=[expand(p) for p in self.support]
        self.disks=[expand(p) for p in self.disks]
        self.rows=[(np.concatenate((a,np.zeros(2,dtype=object))),b) for a,b in self.rows]
        self.numeric_rows=[numeric(a) for a,_ in self.rows]
        self.support.append(eye[previous:previous+2])
        new_point=self.support[-1]

        n,m,c3,_,_,u,_,t,_=setup()
        axis=n[5]
        matrix=2*np.outer(axis,axis)-np.eye(2,dtype=object)
        reflected=matrix@new_point
        for v in n+[m[3],m[5],matrix@m[3],matrix@m[5]]:
            self.linear(np.asarray(v,dtype=object)@new_point,HALF)
        self.linear(np.asarray(m[5],dtype=object)@(eye[6:8]-new_point),0)
        self.linear(np.asarray(m[5],dtype=object)@(eye[8:10]-reflected),0)
        for point in self.witnesses:
            i=len(self.disks)
            self.disks.append(point-new_point)
            for g in ((1,0),(-1,0),(0,1),(0,-1)):
                self.cut(i,g)
        self.witnesses.append(new_point)
        self.disk_num=np.stack([numeric(matrix) for matrix in self.disks])

        sigma=arb("1.5494009336023418793")*PI/180
        for v in (unit(PI/6-sigma),unit(11*PI/6-sigma)):
            for point in self.witnesses:
                self.linear(np.asarray(v,dtype=object)@point,HALF)
                self.linear(np.asarray(v,dtype=object)@(matrix@point),HALF)

        if avoid_ah:
            co,si=unit(arb(delta))
            radial=minus(u,c3)
            u1=plus(c3,minus(scale(co,radial),scale(si,rot(radial))))
            corners=[(x,y) for x in (u1[0],t[0]) for y in (u1[1],t[1])]
            _,outer=polygon(delta,4)
            for _,v,b in outer:
                if all(dot(v,corner)<b for corner in corners):
                    for point in self.witnesses:
                        self.linear(np.asarray(v,dtype=object)@point,b)
                        self.linear(np.asarray(v,dtype=object)@(matrix@point),b)


def root_clearance_objective(program,phi,normal,cap,multiplier):
    n=program.normals
    root3=arb(3).sqrt()
    d0,d1=normal[0]-normal[1]/root3,2*normal[1]/root3
    mu=exact_fraction(multiplier)
    weights=[arb(1),mu-d0,-mu-d1,mu]
    rhs=cap-(d0+d1-mu)/2
    directions=[normal]+list(n[:3])
    for i,w in enumerate(weights):
        if w<0:
            weights[i]=-w
            directions[i]=scale(-1,directions[i])
            rhs-=w
        elif not w>0:
            assert w.is_zero()
    angle=exact_fraction(phi)*PI/180-4*PI/3
    directions=[rotate_angle(v,angle) for v in directions]
    return program.objective(directions,weights),rhs


def spot_phi(program,phi,sign):
    objective=phi_objective(program,phi,sign)
    return program.validate(objective,program.optimize(objective))-arb(3)/2


def spot_edge(program,phi,normal,cap,multiplier):
    obj,threshold=root_clearance_objective(program,phi,normal,cap,Fraction(multiplier))
    return threshold-program.validate(obj,program.optimize(obj))


def certify(program,right,divisions,replay=None):
    left=Fraction(30)
    nodes=[left+(right-left)*i/divisions for i in range(divisions+1)]
    angle_cos=(exact_fraction((right-left)/divisions)*PI/360).cos()
    assert angle_cos>0
    end_obj=phi_objective(program,right,1)
    end_dual=program.optimize(end_obj) if replay is None else replay["right_end_dual"]
    end_bound=program.validate(end_obj,end_dual)-arb(3)/2
    assert end_bound<0, str(end_bound)

    inner,_=polygon(program.delta,4)
    hex_edges={"B1-Q","P-F1","F1-H","E3-D1","D1-G","C3-B1"}
    assert len([edge for edge in inner if edge[0] in hex_edges])==6
    multiplier_options=tuple(map(Fraction,("0","-1","-0.5","1","0.5","-1.5","-2","1.5","2")))
    cell_certs={}
    saved=0
    all_bounds=[]
    for label,normal,cap in inner:
        if label in hex_edges:
            continue
        cache={}

        def dual_bound(j,mu):
            key=f"{j}_{mu}"
            if key not in cache:
                objective,threshold=root_clearance_objective(program,nodes[j],normal,cap,mu)
                dual=program.optimize(objective) if replay is None else replay["cell_certs"][label]["duals"][key]
                cache[key]=(program.validate(objective,dual),threshold,dual)
            return key,cache[key]

        bounds=[]
        chosen=[]
        used=set()
        for j in range(divisions):
            options=multiplier_options if replay is None else (Fraction(replay["cell_certs"][label]["multipliers"][j]),)
            for mu in options:
                key1,(u,cap_shift,_)=dual_bound(j,mu)
                key2,(v,_,_)=dual_bound(j+1,mu)
                support=max(u,v)
                assert support>0
                margin=cap_shift-support/angle_cos
                if margin>0:
                    bounds.append(margin)
                    chosen.append(str(mu))
                    used.update((key1,key2))
                    break
            else:
                raise ArithmeticError(f"clearance unproved for {label} over {nodes[j]}--{nodes[j+1]}")
        result=min(bounds,key=lambda x:x.lower())
        print(label,str(result.lower()),flush=True)
        all_bounds.append((result,label))
        cell_certs[label]={"multipliers":chosen,"duals":{key:cache[key][2] for key in sorted(used)}}
        saved+=len(used)
    min_margin,label=min(all_bounds,key=lambda row:row[0].lower())
    return {"summary":{"delta":program.delta,"angle_root_range_degrees":[str(left),str(right)],
                        "F_right_upper":str(end_bound),"angle_cells":divisions,"arithmetic_precision_bits":ctx.prec,
                        "nonhex_slack_strict_lower":str(min_margin.lower()),"active_edge":label,
                        "tested_nonhex_edges":len(all_bounds),"total_validated_duals":saved+1,
                        "adaptive_cut_rows":len(program.adaptive_cuts),
                        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            "adaptive_cuts":program.adaptive_cuts,"right_end_dual":end_dual,"cell_certs":cell_certs}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("delta")
    parser.add_argument("right_endpoint_degrees")
    parser.add_argument("output",type=Path)
    parser.add_argument("--angle-cells",type=int,default=12)
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
    result=certify(program,Fraction(args.right_endpoint_degrees),args.angle_cells,replay)
    result["summary"]["cut_gain_area_interval"]=geometry["gain_area"]
    args.output.write_text(json.dumps(result,separators=(",",":"))+"\n")
    print(json.dumps(result["summary"],indent=2),flush=True)
