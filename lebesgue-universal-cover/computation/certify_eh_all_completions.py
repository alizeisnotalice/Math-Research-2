#!/usr/bin/env python3
"""Certify individualized placements for all completions of the twelve-point EH witness."""
import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path

from audit_eh_continuous_placement import polygon
from audit_optimized_eh_chord import (
    HALF, ONE, PI, arb, circle_pair, cross, ctx, dot, minus, neg, plus, rot, scale, setup, unit,
)


def norm_squared(point):
    return dot(point,point)


class Lens:
    def __init__(self, point_map):
        self.points=list(point_map.values())
        self.vertices=[]
        for x,y in combinations(self.points,2):
            distance=norm_squared(minus(x,y))
            assert 0<distance<4
            for sign in (-1,1):
                vertex=circle_pair(x,y,sign)
                if not any(norm_squared(minus(vertex,z))>1 for z in self.points):
                    self.vertices.append(vertex)
        self.differences=[]
        for x in self.points:
            row=[minus(x,y) for y in self.points if y is not x]
            self.differences.append([(d,norm_squared(d)) for d in row])

    def upper_support(self, direction):
        # The true direction has norm one even when its Arb enclosure is wide.
        candidates=[dot(direction,v) for v in self.vertices]
        for center,differences in zip(self.points,self.differences):
            if not any(2*dot(direction,d)+d2>0 for d,d2 in differences):
                candidates.append(ONE+dot(direction,center))
        return max(value.upper() for value in candidates)

    def completion_support(self, direction):
        return (ONE-self.upper_support(neg(direction))).union(self.upper_support(direction))


def raw_lens_alignment(point_map):
    # Only the seven exact unit-distance pairs can constrain these normal cones.
    n,m,*_=setup()
    b_to_w=minus(point_map["b"],point_map["w"])
    reflected=neg(minus(point_map["b_Q"],point_map["w_Q"]))
    reflected_m=neg(minus(point_map["w_opp_Q"],point_map["w_Q"]))
    for ray in (b_to_w,m[5]):
        assert ray[0]>0 and ray[1]<0
    for ray in (reflected,reflected_m):
        assert ray[0]>0 and ray[1]>0 and cross(n[1],ray)>0
    return {"b_minus_w":list(map(str,b_to_w)),"m5":list(map(str,m[5])),
            "reflected_b_minus_w":list(map(str,reflected)),"reflected_m5":list(map(str,reflected_m)),
            "no_width_one_directions_in_degrees_open_modulo_180":[0,60]}


def sign_test(lens, degrees):
    phi=degrees*PI/180
    sides=[lens.completion_support(unit(phi+i*PI/3)) for i in range(3)]
    return sides[0]-sides[1]+sides[2]-HALF


def placement(lens, delta, left, right, num_cells, subdivisions):
    n,*_=setup()
    inner,_=polygon(delta,subdivisions)
    hex_edges={"B1-Q","P-F1","F1-H","E3-D1","D1-G","C3-B1"}
    assert len([edge for edge in inner if edge[0] in hex_edges])==6
    nonhex=[]
    root3=arb(3).sqrt()
    for name,normal,offset in inner:
        if name not in hex_edges:
            coeff_0=normal[0]-normal[1]/root3
            coeff_1=2*normal[1]/root3
            nonhex.append((name,normal,offset,coeff_0,coeff_1))
    width=(right-left)/num_cells
    rows=[]
    for i in range(num_cells):
        phi=(left+i*width).union(left+(i+1)*width)*PI/180
        co,si=unit(phi-PI/3)

        def into_body(direction):
            return plus(scale(co,direction),scale(si,rot(direction)))

        a0=HALF-lens.completion_support(into_body(n[0]))
        a1=HALF-lens.completion_support(into_body(n[1]))
        bounds=[]
        for name,normal,offset,d0,d1 in nonhex:
            bound=offset-lens.upper_support(into_body(normal))-d0*a0-d1*a1
            bounds.append((bound,name))
        least,name=min(bounds,key=lambda result:result[0].lower())
        if not all(bound>0 for bound,_ in bounds):
            return {"delta":delta,"passes":False,"first_unproved_cell":i,
                    "slack_lower_bound":str(least.lower()),"edge":name}
        rows.append((least.lower(),name,i))
    margin,label,i=min(rows,key=lambda result:result[0])
    return {"delta":delta,"passes":True,"cells":num_cells,"curve_subdivisions":subdivisions,
            "strict_minimum_nonhex_slack_lower_bound":str(margin),"closest_edge":label,"cell":i}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("witness",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("delta",nargs="+",default=["0.0000155","0.005"])
    parser.add_argument("--num-cells",type=int,default=256)
    parser.add_argument("--curve-subdivisions",type=int,default=8)
    args=parser.parse_args()
    witness=json.loads(args.witness.read_text())
    points={name:tuple(map(arb,p)) for name,p in witness["points"].items()}
    assert all(norm_squared(point)<(arb(3)/5)**2 for point in points.values())
    lens=Lens(points)
    left,right=arb("1.44"),arb("1.47")
    left_sign,right_sign=sign_test(lens,left),sign_test(lens,right)
    assert left_sign<0 and right_sign>0
    results=[placement(lens,d,left,right,args.num_cells,args.curve_subdivisions) for d in args.delta]
    for result in results:
        if result["passes"]:
            certified_slack=arb(result["strict_minimum_nonhex_slack_lower_bound"])
            for eps_text in ("8e-8","2e-8"):
                eps=arb(eps_text)
                if left_sign<-9*eps/2 and right_sign>9*eps/2 and certified_slack>15*eps/2:
                    result["certified_approximate_witness_radius"]=eps_text
                    break
    report={"precision_bits":ctx.prec,"active_lens_corners_with_uncertain_included":len(lens.vertices),
            "raw_lens_alignment":raw_lens_alignment(points),"root_angle_range_degrees":["1.44","1.47"],
            "root_compatibility_at_left_exact_bounds":[str(left_sign.lower()),str(left_sign.upper())],
            "root_compatibility_at_right_exact_bounds":[str(right_sign.lower()),str(right_sign.upper())],
            "rigid_motion":"R_(pi/3-phi), with exact translation centering widths along n0,n1",
            "results":results,
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "input_sha256":hashlib.sha256(args.witness.read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
