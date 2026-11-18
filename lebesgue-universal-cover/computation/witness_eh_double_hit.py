#!/usr/bin/env python3
"""Certify a finite witness for two EH hits at delta = 0.0000155."""
import argparse
from itertools import combinations
import hashlib
import json
from pathlib import Path

from audit_optimized_eh_chord import (
    HALF, ONE, PI, arb, circle_pair, cross, ctx, dot, minus, neg, plus, rot, scale, setup, trial, unit,
)


def squared(p):
    return dot(p, p)


def witness(delta):
    # Check the arc order, orientation, and all separate relocation signs.
    trial(delta, tight=True)
    n, m, c3, b3, pal, u, _, t, j = setup()
    c, s = unit(arb(delta))
    radial = minus(u, c3)
    u1 = plus(c3, minus(scale(c, radial), scale(s, rot(radial))))
    direction = minus(t, u1)
    k = direction[0]/(-direction[1])
    gamma = HALF-arb(3).sqrt()*k/2
    height = (HALF-dot(pal, u1))/gamma
    H, a = HALF-u1[0], c3[1]-u1[1]
    cb_limit = (ONE-(H-k*height)**2).sqrt()-a-height
    r = cb_limit-arb("1e-7")
    b = plus(c3, (arb(0),r))
    z_chord = plus(u1, (k*height,-height))
    z = plus(z_chord, (arb("-1e-10"),arb("2e-10")))
    axis = n[5]

    def reflection(x):
        return minus(scale(2*dot(axis,x),axis),x)

    eps = arb("1e-8")
    normal_term = -m[5][1]*r-eps
    assert r>0 and 0<normal_term<1
    tangent_term = (2*normal_term-normal_term**2).sqrt()
    w = minus(plus(b,scale(normal_term-1,m[5])),scale(tangent_term,rot(m[5])))
    points = {"z":z,"b":b,"b_opp":minus(b,n[0]),"w":w,"w_opp":plus(w,m[5])}
    points.update({name+"_Q":reflection(x) for name,x in list(points.items())})
    points.update({"q":scale(HALF,axis),"q_opp":scale(-HALF,axis)})

    # P ∩ QP. Six equalities follow exactly from the construction and Qn_j permutation.
    exact_normals = {"b":0,"b_opp":3,"b_Q":4,"b_opp_Q":1,"q":5,"q_opp":2}
    normals = n+[m[3],m[5],reflection(m[3]),reflection(m[5])]
    planes = [(HALF-dot(v,x),name,index) for name,x in points.items()
              for index,v in enumerate(normals) if exact_normals.get(name)!=index]
    assert all(gap>0 for gap,_,_ in planes)

    # The algebraic equalities have unit difference n0,n5,m5, or use
    # (1-normal_term)^2+tangent_term^2=1, and also their reflected counterparts.
    equal_pairs = {frozenset(pair) for pair in (("b","b_opp"),("b","w"),
                                                ("w","w_opp"),("q","q_opp"))}
    equal_pairs.update({frozenset(name+"_Q" for name in pair) for pair in list(equal_pairs)
                        if "q" not in pair})
    distances = [(ONE-squared(minus(x,y)),name,other)
                 for (name,x),(other,y) in combinations(points.items(),2)
                 if frozenset((name,other)) not in equal_pairs]
    assert all(gap>0 for gap,_,_ in distances)

    # Interior of the EH lens cap and its exact defining ray orientations.
    ear = {
        "distance_to_C3_squared_slack": ONE-squared(minus(z,c3)),
        "distance_to_B3_squared_slack": ONE-squared(minus(z,b3)),
        "chord_right_side_slack": -cross(direction,minus(z,u1)),
        "U1_inside_B3_disk_slack": ONE-squared(minus(u1,b3)),
        "T_inside_C3_disk_slack": ONE-squared(minus(t,c3)),
        "C3_disk_entry_from_U1": -dot(direction,minus(u1,c3)),
        "B3_disk_entry_from_T": dot(direction,minus(t,b3)),
        "J_on_chord_right_side": -cross(direction,minus(j,u1)),
        "other_circle_crossing_on_chord_left_side":cross(direction,minus(circle_pair(b3,c3,1),u1)),
    }
    assert all(value>0 for value in ear.values())

    sigma=arb("1.5494009336023418793")*PI/180
    flags = {
        "A_plus": dot(m[1],points["w_Q"])-HALF,
        "A_minus": dot(unit(PI/2-sigma),points["w_Q"])-HALF,
        "F_minus": dot(unit(5*PI/6-sigma),w)-HALF,
        "opposite_m3": dot(neg(m[3]),points["z_Q"])-HALF,
        "opposite_m5": dot(neg(m[5]),w)-HALF,
    }
    assert all(value>0 for value in flags.values())

    excluded_placements = {}
    for i, (co, si) in enumerate(n):
        for reflect in (False, True):
            if (i,reflect) in ((0,False),(4,True)):
                continue  # Identity and Q2 fail by the interior ear hit.
            cuts = []
            for name,(x,y) in points.items():
                yr = -y if reflect else y
                image = (co*x-si*yr,si*x+co*yr)
                for index in (3,5):
                    cuts.append((dot(m[index],image)-HALF,name,index))
            gap,name,index = max(cuts,key=lambda row: row[0].lower())
            assert gap>0
            excluded_placements[f"rotation_{i}_reflection_{int(reflect)}"] = [str(gap),name,index]

    plane_gap,*plane_at=min(planes,key=lambda row: row[0].lower())
    distance_gap,*distance_at=min(distances,key=lambda row: row[0].lower())
    return {
        "delta":delta, "precision":ctx.prec, "points":{k:list(map(str,x)) for k,x in points.items()},
        "other_halfplane_minimum_slack":[str(plane_gap),*plane_at],
        "other_pair_minimum_squared_distance_slack":[str(distance_gap),*distance_at],
        "ear":{k:str(v) for k,v in ear.items()}, "flags":{k:str(v) for k,v in flags.items()},
        "other_dihedral_placements":excluded_placements,
        "remaining_halfplanes_at_equality":exact_normals,
        "remaining_pairs_at_equality":sorted(sorted(pair) for pair in equal_pairs),
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    result=witness("0.0000155")
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
