#!/usr/bin/env python3
"""Bound the placement angles of the twelve-point witness over a full turn."""
import argparse
from itertools import combinations, groupby
import json
from pathlib import Path

import numpy as np

from audit_eh_continuous_placement import polygon, rigid_check
from audit_optimized_eh_chord import PI, arb, cross, dot, norm, plus, rot, scale, unit


def coordinate_midpoints(points):
    return np.array([[float(c.mid()) for c in point] for point in points])


def linear_data(edges):
    normals=coordinate_midpoints([v for _,v,_ in edges])
    offsets=np.array([float(h.mid()) for _,_,h in edges])
    return normals,offsets


def translation_candidates(normals):
    augmented=np.column_stack((normals,np.ones(len(normals))))
    triples=np.array(list(combinations(range(len(normals)),3)))
    matrices=augmented[triples]
    stable=np.abs(np.linalg.det(matrices))>1e-9
    return augmented,triples[stable],np.linalg.inv(matrices[stable])


def exclusion_candidates(normals):
    triples=np.array(list(combinations(range(len(normals)),3)))
    selected=normals[triples]
    a,b,c=selected[:,0,:],selected[:,1,:],selected[:,2,:]
    det=lambda x,y: x[:,0]*y[:,1]-x[:,1]*y[:,0]
    weights=np.column_stack((det(b,c),det(c,a),det(a,b)))
    positive=(np.min(weights,axis=1)>1e-12)|(np.max(weights,axis=1)<-1e-12)
    weights=weights[positive]
    return triples[positive],weights/np.sum(weights,axis=1)[:,None]


def classify(delta, point_map, count, subdivisions, refine):
    inner,outer=polygon(delta,subdivisions)
    iv,ih=linear_data(inner)
    ov,oh=linear_data(outer)
    im,it,im_inv=translation_candidates(iv)
    ot,ow=exclusion_candidates(ov)
    points=list(point_map.values())
    xp=coordinate_midpoints(points)
    radius=arb(3)/5
    assert all(norm(point)<radius for point in points)
    denominator=count*2**refine
    bins=[]
    proof_rows=[]

    def visit(start,after):
        cell_width=after-start
        angle=PI*(start+after)/denominator
        half_width=PI*cell_width/denominator
        allowance=radius*half_width
        pad=float(allowance.upper())
        a=float(angle.mid())
        co,si=np.cos(a),np.sin(a)
        rotated=xp@np.array(((co,si),(-si,co)))
        od=rotated@ov.T
        best=od.argmax(axis=0)
        outer_remaining=oh-od[best,np.arange(len(ov))]
        estimates=np.sum(ow*outer_remaining[ot],axis=1)
        k=estimates.argmin()
        if estimates[k]<0:
            selected=[outer[i] for i in ot[k]]
            va,vb,vc=(v for _,v,_ in selected)
            weights=[cross(vb,vc),cross(vc,va),cross(va,vb)]
            if all(w<0 for w in weights):
                weights=[-w for w in weights]
            co,si=unit(angle)
            terms=[]
            slopes=[]
            for i,(_,normal,offset) in zip(ot[k],selected):
                point=points[best[i]]
                rotated_point=plus(scale(co,point),scale(si,rot(point)))
                terms.append(offset-dot(normal,rotated_point))
                slopes.append(dot(normal,rot(rotated_point)))
            remainder=abs(sum(w*s for w,s in zip(weights,slopes)))*half_width
            remainder+=sum(weights)*radius*half_width**2/2
            if all(w>0 for w in weights) and sum(w*b for w,b in zip(weights,terms))+remainder<0:
                bins.append([start,after,"excluded"])
                proof_rows.append([start,after,"excluded",[label for label,_,_ in selected],
                                   [list(point_map)[best[i]] for i in ot[k]]])
                return

        rh=ih-(rotated@iv.T).max(axis=0)
        candidates=np.einsum('bij,bj->bi',im_inv,rh[it])
        feasible=np.max(candidates@im.T-rh,axis=1)<1e-10
        present=candidates[feasible]
        if len(present):
            t=present[np.argmax(present[:,2])]
            if t[2]>pad:
                tx,ty=(format(value,'.14g') for value in t[:2])
                co,si=unit(angle)
                shifted=[plus((arb(tx),arb(ty)),plus(scale(co,p),scale(si,rot(p)))) for p in points]
                if all(offset-dot(normal,x)>allowance for x in shifted for _,normal,offset in inner):
                    bins.append([start,after,"feasible"])
                    proof_rows.append([start,after,"feasible",[tx,ty]])
                    return
        if cell_width>1:
            middle=(start+after)//2
            visit(start,middle)
            visit(middle,after)
        else:
            bins.append([start,after,"undecided"])

    stride=2**refine
    for start in range(0,denominator,stride):
        visit(start,start+stride)

    bands=[]
    counts={"excluded":0,"feasible":0,"undecided":0}
    for status,entries in groupby(bins,key=lambda row:row[2]):
        group=list(entries)
        first,after=group[0][0],group[-1][1]
        bands.append([first,after,status])
        counts[status]+=after-first
    spot=rigid_check(inner,point_map,"4.489","-0.0156","0.0294")
    assert spot["contained"]
    assert sum(counts.values())==denominator
    return {"delta":delta,"angle_cell_denominator_in_turns":denominator,"curve_subdivisions":subdivisions,
            "counts":counts,
            "bands_in_turn_numerators":bands,"strict_placement":spot,"proof_rows":proof_rows}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("witness",type=Path)
    parser.add_argument("delta")
    parser.add_argument("count",type=int)
    parser.add_argument("output",type=Path)
    parser.add_argument("--curve-subdivisions",type=int,default=4)
    parser.add_argument("--refine",type=int,default=0)
    args=parser.parse_args()
    points={name:tuple(map(arb,p)) for name,p in json.loads(args.witness.read_text())["points"].items()}
    report=classify(args.delta,points,args.count,args.curve_subdivisions,args.refine)
    args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k!="proof_rows"},indent=2))
