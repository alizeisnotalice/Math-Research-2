#!/usr/bin/env python3
"""Certify the coupled double-hit support bounds by supported-disk linear relaxations."""
import argparse
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from audit_eh_continuous_placement import polygon
from audit_optimized_eh_chord import (
    HALF, ONE, PI, arb, ctx, dot, minus, plus, rot, scale, setup, trial, unit,
)


def midpoint(x):
    return float(x.mid()) if isinstance(x,arb) else float(x)


def numeric(array):
    return np.vectorize(midpoint,otypes=[float])(array)


class Program:
    def __init__(self,delta):
        n,m,c3,b3,p,u,_,t,_=setup()
        co,si=unit(arb(delta))
        rad=minus(u,c3)
        u1=plus(c3,minus(scale(co,rad),scale(si,rot(rad))))
        chord=minus(t,u1)
        axis=n[5]
        r3=arb(3).sqrt()
        self.normals=n
        self.delta=delta
        self.dimension=17
        coordinates=np.eye(self.dimension+1,dtype=object)
        constant=lambda x: np.asarray(x,dtype=object)[:,None]*coordinates[-1]
        vec=lambda i: coordinates[i:i+2]
        matrix=2*np.outer(axis,axis)-np.eye(2,dtype=object)
        reflect=lambda x: matrix@x
        z,zp=vec(0),vec(2)
        b,bp=(constant(c3)+np.outer((0,1),coordinates[i]) for i in (4,5))
        w,wp=vec(6),vec(8)
        q=constant(scale(HALF,axis))+np.outer(rot(axis),coordinates[10])
        self.support=[vec(i) for i in (11,13,15)]
        self.rows=[]
        self.numeric_rows=[]
        self.numeric_rhs=[]
        self.adaptive_cuts=[]
        self.disks=[]
        self.witnesses=[z,b,b-constant(n[0]),w,w+constant(m[5]),q,q-constant(axis)]
        self.witnesses.extend(reflect(y) for y in (zp,bp,bp-constant(n[0]),wp,wp+constant(m[5])))
        self.witnesses.extend(self.support)
        for point in self.witnesses:
            for normal in n+[m[3],m[5],matrix@m[3],matrix@m[5]]:
                self.linear(np.asarray(normal,dtype=object)@point,HALF)
            self.linear(np.asarray(m[5],dtype=object)@(w-point),0)
            self.linear(np.asarray(m[5],dtype=object)@(wp-reflect(point)),0)

        # Both contacts are on the C3 side of the CB line. A' is the Q image of -m5.
        self.linear(-coordinates[4],0)
        self.linear(-coordinates[5],0)
        self.linear(np.asarray(m[5],dtype=object)@w,-HALF)
        self.linear(np.asarray(m[5],dtype=object)@wp,-HALF)
        for hit in (z,zp):
            self.linear(np.asarray(rot(chord),dtype=object)@(hit-constant(u1)),0)
            self.linear(hit[1]-u1[1]*coordinates[-1],0)
            self.disks.extend(hit-constant(center) for center in (c3,b3))
        self.disks.extend(a-b for a,b in combinations(self.witnesses,2))
        self.disk_num=np.stack([numeric(matrix) for matrix in self.disks])
        for k in range(len(self.disks)):
            for g in ((1,0),(-1,0),(0,1),(0,-1)):
                self.cut(k,g)
        self.initial_rows=len(self.rows)

    def linear(self,coeffs,bound):
        self.rows.append((coeffs[:-1].copy(),bound-coeffs[-1]))
        self.numeric_rows.append(numeric(self.rows[-1][0]))
        self.numeric_rhs.append(midpoint(self.rows[-1][1]))

    def cut(self,i,g,record=False):
        if record:
            self.adaptive_cuts.append([int(i),[float(x) for x in g]])
        if any(isinstance(x,float) for x in g):
            g=tuple(arb(x) for x in g)
        self.linear(np.asarray(g,dtype=object)@self.disks[i],arb(dot(g,g)).sqrt())

    def objective(self,directions,weights):
        return sum(np.asarray(u,dtype=object)@x*w for u,x,w in zip(directions,self.support,weights))

    def optimize(self,objective,tolerance=2e-10):
        from scipy.optimize import linprog

        c=numeric(objective[:-1])
        for _ in range(45):
            solution=linprog(-c,A_ub=self.numeric_rows,b_ub=self.numeric_rhs,bounds=(-1,1),
                             method="highs",options={"primal_feasibility_tolerance":1e-10,
                                                     "dual_feasibility_tolerance":1e-10})
            if not solution.success:
                raise ValueError(solution.message)
            disk_offsets=self.disk_num@np.append(solution.x,1.)
            lengths=np.linalg.norm(disk_offsets,axis=1)
            if max(lengths)<=1+tolerance:
                break
            indices=np.argsort(lengths)[-24:]
            for i in indices[lengths[indices]>1+tolerance]:
                self.cut(i,tuple(disk_offsets[i]/lengths[i]),record=True)

        dual=[(int(i),float(-dual)) for i,dual in enumerate(solution.ineqlin.marginals) if dual<0]
        return dual

    def validate(self,objective,dual):
        # Any nonnegative inequality weights give a rigorous linear objective upper bound.
        residual=np.asarray(objective[:-1],dtype=object).copy()
        upper=arb(objective[-1])
        for i,weight in dual:
            assert 0<=i<len(self.rows) and weight>0
            mul=arb(weight)
            ca,cb=self.rows[i]
            residual-=ca*mul
            upper+=cb*mul
        upper+=sum(abs(value) for value in residual)  # every physical coordinate has abs <= 1
        return upper.upper()


def rotate_angle(point,angle):
    c,s=unit(angle)
    return plus(scale(c,point),scale(s,rot(point)))


def exact_fraction(f):
    return arb(f.numerator)/f.denominator


def phi_objective(program,degrees,sign):
    directions=[unit(exact_fraction(degrees)*PI/180+j*PI/3) for j in range(3)]
    directions=[scale(f*sign,n) for n,f in zip(directions,(1,-1,1))]
    return program.objective(directions,(1,1,1))


def edge_objective(program,degrees,normal,cap):
    n=program.normals
    root3=arb(3).sqrt()
    factors=(normal[0]-normal[1]/root3,2*normal[1]/root3)
    sv=[normal]
    weights=[arb(1)]
    for d,base in zip(factors,n):
        if d>0:
            sv.append(scale(-1,base))
            weights.append(d)
        elif d<0:
            sv.append(base)
            weights.append(-d)
        else:
            assert d.is_zero()
            sv.append(base)
            weights.append(abs(d))
    premove=exact_fraction(degrees)*PI/180-PI/3
    directions=[rotate_angle(v,premove) for v in sv]
    return program.objective(directions,weights),cap+(weights[1]+weights[2])/2


def certify(model,saved_duals=None):
    left,right=Fraction("0.65"),Fraction("2.0")
    divisions=14
    nodes=[left+(right-left)*i/divisions for i in range(divisions+1)]
    duals={}

    def compute(key,objective):
        dual=model.optimize(objective) if saved_duals is None else saved_duals[key]
        duals[key]=dual
        return model.validate(objective,dual)

    left_upper=compute("F_left",phi_objective(model,left,1))-arb(3)/2
    right_neg_upper=compute("negative_F_right",phi_objective(model,right,-1))-arb(3)/2
    assert left_upper<0 and right_neg_upper<0

    inner,_=polygon(model.delta,4)
    hex_edges={"B1-Q","P-F1","F1-H","E3-D1","D1-G","C3-B1"}
    assert len([edge for edge in inner if edge[0] in hex_edges])==6
    gain=(exact_fraction((right-left)/divisions)*PI/360).cos()
    assert gain>0
    bounds=[]
    for label,normal,cap in inner:
        if label in hex_edges:
            continue
        values=[]
        for j,node in enumerate(nodes):
            objective,threshold=edge_objective(model,node,normal,cap)
            upper=compute(label+"_"+str(j),objective)
            values.append((upper,threshold))
        candidates=[threshold-max(u,v)/gain for (u,threshold),(v,_) in zip(values,values[1:])]
        margin=min(candidates,key=lambda x:x.lower())
        assert margin>0,(label,str(margin))
        bounds.append((margin,label))
        print(label,str(margin.lower()),flush=True)
    if saved_duals is not None:
        assert set(saved_duals)==set(duals)
    margin,label=min(bounds,key=lambda row:row[0].lower())
    result={"delta":model.delta,"arithmetic_precision_bits":ctx.prec,
            "compatibility_root_range_degrees":[str(left),str(right)],
            "F_left_upper":str(left_upper), "F_right_lower":str(-right_neg_upper),
            "angle_interpolation_cells":divisions,"inner_polygon_curve_subdivisions":4,
            "bounded_nonhex_edges":len(bounds),"required_duals":len(duals),
            "nonhex_slack_strict_lower":str(margin.lower()),"active_edge":label,
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    return result,duals


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("delta")
    parser.add_argument("output",type=Path)
    parser.add_argument("--replay",type=Path)
    args=parser.parse_args()
    geometry=trial(args.delta,tight=True)
    model=Program(args.delta)
    saved_duals=None
    if args.replay:
        saved=json.loads(args.replay.read_text())
        assert saved["summary"]["delta"]==args.delta
        for i,g in saved["adaptive_cuts"]:
            model.cut(i,g,record=True)
        saved_duals=saved["duals"]
    result,duals=certify(model,saved_duals)
    result["adaptive_cut_rows"]=len(model.adaptive_cuts)
    result["cut_gain_area_interval"]=geometry["gain_area"]
    print(json.dumps(result,indent=2),flush=True)
    args.output.write_text(json.dumps({"summary":result,"adaptive_cuts":model.adaptive_cuts,"duals":duals},
                                     separators=(",",":"))+"\n")
