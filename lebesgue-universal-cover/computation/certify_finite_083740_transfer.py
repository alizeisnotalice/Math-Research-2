"""Replay the shrink transfer at the R161 threshold and construct diameter-one rational groups."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import platform
from zipfile import ZipFile

import flint
from flint import arb, ctx
import numpy as np

from certify_finite_point_leaf import point_tables
from certify_nine import PI_UPPER, ar, hull
from certify_polygon_083731 import (
    HERE, LAMBDA, ORDERS, diameter_squared, lower_radius, sampled_group,
    verified_table, verify_core_table,
)
from certify_tomar_six_high import SourceTree
from select_tomar_083731_scores import AREA, radius_rows
from tangent_polygons_finite_cover import tangent_points

TARGET = F(".83740")
GOAL = F(int(TARGET*AREA)+1, AREA)
EXTRA_TARGET = F(".837403")
RHO = F(".000000555")
FOUR = {(1205, pos) for pos in (139716,139730,151617,152529)}
PINNED = {
    "finite-point-lift-083740.json": "267c8e2dd03a5f914908ee31f04b68c3762cb92c109be5a3f2edb22d734dde24",
    "tomar-083740-score-layers-audit.json": "87d4f754768c8723461c20cfa88f95a5a4deda06dea5b58aa659fd86e37b8b0b",
    "joint-contraction-audit.json": "6f86deeb0a70431f3e299138ff3483139aa9a2b008b59945a7a2d9d4af37206e",
}
FOUR_KIT_SHA = "0d0810dc095a29ace19741bc1490adf19c6b7a00550ec90ff41d28e81104f8c6"
FOUR_MAC_SHA = "b43e2de7409ae48d696bc3fff94759827b92bb66f86d9239974c8a8c51af79d8"


def checksum(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path, expected=None):
    if expected is not None:
        assert checksum(path) == expected
    return json.loads(path.read_bytes())


def named(name, expected=None):
    return read(HERE/name, PINNED.get(name, expected))


def audit_cohorts(a, tree):
    terminal = named("finite-point-lift-083740.json")
    stage = named("tomar-083740-stage-ledger.json",terminal["outer_audit_sha256"])
    audit = named("tomar-083740-score-layers-audit.json")
    seed = named("tomar-083740-seed.json",stage["seed_receipt_sha"])
    retained = named("tomar-083740-retained.json",terminal["retained_sha256"])
    assert F(terminal["target"]) == F(stage["target"]) == F(audit["target"]) == F(seed["target"]) == TARGET
    assert stage["integer_goal"] == seed["target_integer"] == GOAL*AREA
    assert stage["complete"] and stage["parts"] == len(seed["rows"]) == 1024
    assert audit["retained"] == retained
    assert audit["stage_sha256"] == terminal["outer_audit_sha256"]
    assert audit["first_seed_sha256"] == stage["seed_receipt_sha"]
    assert checksum(a.tree/"SHA256SUMS") == seed["tree_manifest_sha256"]
    assert F(seed["radial_target"]) == F(stage["radial_target"]) == F(".83741")
    assert len(stage["audits"]) == 4
    remain = seed["remaining"]
    for row in stage["audits"]:
        assert remain == row["tested"] == row["accepted"]+row["remaining"]
        remain = row["remaining"]
    all_keys = {(int(p),pos) for p,ids in retained.items() for pos in ids}
    assert len(all_keys) == remain == 1940
    rows = terminal["rows"]
    assert len(rows) == 1936 and [tuple(r[:2]) for r in rows] == sorted(all_keys-FOUR)
    for _,_,n,calls,old,new,depth,tried,joint in rows:
        assert n in (32,64) and 0<=joint<=new<=tried and depth>=0
        assert calls==2*(old+new)-1
    volume = {k:sum(r[j] for r in rows) for k,j in
              (("nodes",3),("five_terminals",4),("point_terminals",5),("joint_terminals",8))}
    volume["point_used"] = sum(r[5]>0 for r in rows)
    assert volume == {"nodes":610999790,"five_terminals":270747453,
                      "point_terminals":34753410,"joint_terminals":34753410,"point_used":1136}

    if a.ordinary_logs:
        ordinary = {tuple(row[:2]):row for row in rows}
        reproduced = {}
        raw = a.ordinary_logs
        logs = [json.loads(path.read_bytes()) for path in raw.glob("weighted-*.json")]
        logs += [json.loads(line) for line in (raw/"stream-extra.jsonl").read_text().splitlines()]
        for row in logs:
            key = (row["part"],row["pos"])
            if key not in ordinary:
                assert key in FOUR
                continue
            assert key not in reproduced
            _,_,n,calls,old,new,depth,tried,joint = ordinary[key]
            status = row["stat"]
            assert row.get("table",row.get("level")) == n
            assert status[:7] == [1,calls,old,new,depth,0,tried] and status[29] == joint
            reproduced[key] = status[:7]+[status[29]]
        assert reproduced.keys() == ordinary.keys()
        volume["ordinary_raw_records_checked"] = len(reproduced)

    padded = read(a.four_macbook,FOUR_MAC_SHA)
    assert checksum(a.four_kit) == FOUR_KIT_SHA
    assert padded["exact_branch_audit"] and padded["exit"]==0 and padded["complete_four"]
    assert F(padded["target"]) == EXTRA_TARGET
    assert {(padded["part"],int(k)) for k in padded["leaves"]} == FOUR
    assert sum(row["visits"] for row in padded["leaves"].values()) == padded["nodes"] == 7209
    with ZipFile(a.four_kit) as kit:
        for key, row in padded["leaves"].items():
            stem = f"part1205-pos{key}"
            packed = kit.read("expected/"+stem+".json")
            payload = kit.read("expected/"+stem+".bin")
            assert sha256(packed).hexdigest()==row["record_sha256"]
            assert sha256(payload).hexdigest()==row["trace_sha256"]
            case = json.loads(packed)
            assert F(case["threshold"])==EXTRA_TARGET and (case["part"],case["pos"])==(1205,int(key))
            assert case["upstream"]==padded["upstream_revision"] and not case["result"]["remaining"]
            assert case["original_header_sha256"]==tree.pins["part1205.json"]
            assert case["original_tree_sha256"]==tree.pins["part1205.npy"]
            box=tree.leaf(1205,int(key))
            assert [z*2**44 for side in zip(*box,strict=True) for z in side]==case["root"]
            if a.local_four:
                assert (a.local_four/(stem+".json")).read_bytes()==packed
                assert (a.local_four/(stem+".bin")).read_bytes()==payload
    radial = named("five-body-radial-affine.json",seed["provenance"]["radial_certificate"])
    assert radius_rows(radial,F(seed["radial_target"])) == seed["radii"]
    assert all(F(r["integer_area_lower"],AREA)>GOAL for r in seed["radii"])
    index = named("tomar-multi-target-after-083735.json",seed["provenance"]["whole_leaf_index"])
    assert index["area_denominator"] == AREA and index["canonical_leaves"] == 6710314
    parent_index = named("tomar-multi-target-scores.json")
    # The only extra overall five-family floor in this index is below the gate.
    assert parent_index["prior_global_integer"] < GOAL*AREA
    layer_totals = [{k:layer[k] for k in ("tested","accepted","remaining","calls")}
                    for layer in stage["audits"]]
    return terminal, radial, {
        "ordinary":volume, "outer_layers":layer_totals, "ordinary_leaves":len(rows),
        "outer_retained":len(all_keys), "dynamic_leaves":len(padded["leaves"]),
        "four_nodes":padded["nodes"], "four_local_exact_match":a.local_four is not None,
        "source_hashes":{**PINNED,"retained":terminal["retained_sha256"],
                         "outer_ledger":terminal["outer_audit_sha256"],
                         "outer_seed":stage["seed_receipt_sha"],"four_kit":FOUR_KIT_SHA,
                         "four_cross_host":FOUR_MAC_SHA}}


def audit_geometry(a, tree, terminal, radial):
    pi = arb.pi()
    assert pi < ar(PI_UPPER)
    cloud = named("finite-point-leaf-input.json",terminal["point_sha256"])
    bank = named("finite-point-global-input.json",terminal["bank_sha256"])
    base = named("tomar-six-high-input.json",terminal["base_sha256"])
    assert checksum(HERE/"tomar_six_083730_guard.cpp") == terminal["weighted_sha256"]
    d2, ptr, triples = point_tables({**cloud,"cells":bank["cells"]})
    assert len(ptr)==385 and len(triples)//8==2493 and d2==100000000008865717121
    tables = {}
    for n in (32,64):
        meta = base["tables"][str(n)]
        cnt,pts,stamp=verified_table(getattr(a,f"table{n}"),tree,n,meta)
        verify_core_table(tree,cnt,pts,pi)
        tables[str(n)] = {"npz_sha256":stamp,"canonical_arrays_sha256":meta["canonical_arrays_sha256"],
                          "shape":list(pts.shape)}
        if n==32:
            assert stamp == radial["source_table_sha256"]
            path = HERE/"five-body-radial-material.npz"
            assert checksum(path) == radial["table_sha"]
            assert checksum(HERE/"certify_radial_affine_fans.py") == radial["source_sha256"]
            with np.load(path,allow_pickle=False) as fixed:
                assert np.array_equal(fixed["counts"],cnt) and np.array_equal(fixed["points"],pts[0])

    r = F(1,2)
    r_small = LAMBDA*r
    R = F(".703")
    rad_a = ar(r_small)
    R_a = ar(R)
    spike = (pi-(rad_a/R_a).acos())*rad_a**2+rad_a*(R_a**2-rad_a**2).sqrt()
    assert spike > ar(GOAL) and R**2 < 2*r_small**2
    root_bounds = []
    for b,n in enumerate(ORDERS[1:]):
        j=3*b
        distance=min(F(tree.root["hi"][j]),F(tree.root["hi"][j+1]),
                     -F(tree.root["lo"][j]),-F(tree.root["lo"][j+1]))
        Rn = ar(LAMBDA)*lower_radius(n,pi)
        val=pi*(rad_a**2+Rn**2)/2+ar(distance)*(rad_a+Rn)
        assert val > ar(F(".859")) > ar(GOAL)
        root_bounds.append(str(val))

    # Replicate the full-depth R152 comparison using exact fractions.
    a0 = F(".42264973")
    assert all(lower_radius(n,pi) > ar(a0) for n in ORDERS)
    point_radius = F(".528619")
    eps = 1-LAMBDA
    point_loss = RHO*point_radius
    ratio = point_loss/eps
    for x,y in cloud["points"]:
        den=cloud["denominator"]
        assert (F(x,den)+F(".01007"))**2+F(y,den)**2 < point_radius**2
    max_depth = F(".204")
    assert R + point_loss-r < max_depth
    v = (1+LAMBDA)*r*r+ratio**2*eps
    assert (1+LAMBDA)-(2*ratio*r+v)/(4*r*r)>0 and ratio<a0
    upper = (1+LAMBDA)*max_depth + max(F(0),(v-2*ratio*(r+max_depth))/(2*r+max_depth))
    margin = 4*(a0-ratio)*LAMBDA**2-upper
    assert margin > F(".0067")
    known = named("joint-contraction-audit.json")
    assert F(known["rho"]) == RHO and F(known["comparison_margin"]) == margin

    # For the four dynamic certificates, use geometric stability rather than
    # assuming their own material tables are in the contracted bodies.
    max_radius = F(3,5)
    assert all(1/(2*(pi/(2*n)).cos())<ar(max_radius) for n in ORDERS[1:])
    assert r < max_radius
    delta = eps*max_radius
    steiner = PI_UPPER*(2*R*delta+delta**2)
    strict = EXTRA_TARGET-GOAL-steiner
    assert strict > F(".00000044")
    return cloud, dict(tables=tables, point_bank=[384,2493], point_diameter_squared_num=d2,
                       cap_comparison_margin=str(margin), radius_area_lower=str(spike),
                       root_exterior_area_lowers=root_bounds, max_old_cap_depth=str(max_depth),
                       dynamic_loss_upper=str(steiner), dynamic_margin=str(strict))


def scale_sixth(cloud, f):
    s=f*(1-RHO)
    return {"points":[[x*s.numerator,y*s.numerator] for x,y in cloud["points"]],
            "denominator":cloud["denominator"]*s.denominator}


def output_groups(path, groups, f):
    displayed = [diameter_squared(g) for g in groups]
    assert all(d<=1 for d in displayed)
    floor = GOAL*f*f
    assert floor > F(".83740093")
    data = (json.dumps(groups,separators=(",",":"))+"\n").encode()
    path.write_bytes(data)
    sizes = [len(g["points"]) for g in groups]
    return dict(points_sha256=sha256(data).hexdigest(),sizes=sizes,total=sum(sizes),
                common_factor=str(f),enclosed_core_scale=str(f*LAMBDA),
                rational_area_floor=str(floor),area_floor_display=float(floor),
                squared_diameters=list(map(str,displayed)))


def certified_core_scale(n, pts, den):
    # Independently verify each rational polygon edge against support of the
    # whole centered curved body. No grid-of-directions approximation enters.
    poly=hull(pts)
    assert len(poly)>=3
    pi=arb.pi()
    corners=[]
    if n:
        rad=1/(2*(pi/(2*n)).cos())
        for j in range(n):
            sin,cos=(2*pi*j/n).sin_cos()
            corners.append((rad*cos,rad*sin))
    lower=arb(2)
    for a,b in zip(poly,poly[1:]+poly[:1],strict=True):
        u,v=b[1]-a[1],a[0]-b[0]
        length=arb(u*u+v*v).sqrt()
        wall=u*a[0]+v*a[1]
        assert wall>0 and length>0
        if not n:
            support=length/2
        else:
            support=arb(0)
            for x,y in corners:
                dot=u*x+v*y
                support=max(support,dot.upper())
                # This radial unit-arc point is admissible exactly when the
                # opposite direction is within pi/(2n) of its center vertex.
                if not dot > -length/2:
                    support=max(support,(dot+length).upper())
        assert support>0
        lower=min(lower,(arb(wall)/(support*den)).lower())
    grid=10**18
    num=int((lower*grid).lower().floor().unique_fmpz())
    scale=F(num,grid)
    assert 0<ar(scale)<lower
    return scale,str(lower),len(poly)


def tangent_variant(m_array, f_max, cloud):
    groups=[]
    fits=[]
    f=f_max
    for n,m in zip(ORDERS,m_array,strict=True):
        pts=tangent_points(n,m)
        assert len(pts)==(n*(m+1) if n else m)
        D=max((x-a)**2+(y-b)**2 for (x,y),(a,b) in combinations(pts,2))
        den=isqrt(D-1)+1
        assert (den-1)**2 < D <= den**2
        kappa,arb_bound,edge_count=certified_core_scale(n,pts,den)
        f=min(f,kappa/LAMBDA)
        fits.append(dict(core_scale=str(kappa),arb_scale_lower=arb_bound,polygon_edges=edge_count))
        groups.append(dict(body=n,denominator=den,points=pts))
    assert 0<LAMBDA*f<1
    assert all(LAMBDA*f <= F(row["core_scale"]) for row in fits)
    groups.append(scale_sixth(cloud,f))
    return groups,f,fits


def point_count_limit():
    # The published vertex-count bounds apply to every polygon in this
    # particular containment-and-scaling transfer at the prescribed target.
    pi=arb.pi()
    target=F(".83740093")
    dilation=1/(ar(LAMBDA)*ar(target/GOAL).sqrt())
    disk_gaps=[(pi/(2*m)).cos()/(pi/m).cos()-dilation for m in (3005,3006)]
    assert disk_gaps[0]>0>disk_gaps[1]
    count_lower=[3006]
    arc_intervals=[]
    for n in ORDERS[1:]:
        required=pi/(2*n*(1/dilation).acos())
        num=int(required.lower().floor().unique_fmpz())
        assert num<required<num+1
        arc_intervals.append(str(required))
        count_lower.append(n*num)
    count_lower.append(30)
    assert sum(count_lower)==9962
    return dict(target=str(target),maximum_dilation=str(dilation),disk_gaps=list(map(str,disk_gaps)),
                arc_intervals=arc_intervals,per_group=count_lower,total=sum(count_lower))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("tree","table32","table64","four_kit","four_macbook","out"):
        p.add_argument(name,type=Path)
    p.add_argument("--local-four",type=Path)
    p.add_argument("--ordinary-logs",type=Path)
    a=p.parse_args()
    assert __debug__
    ctx.prec=160
    tree=SourceTree(a.tree)
    terminal,radial,counts=audit_cohorts(a,tree)
    cloud,geometry=audit_geometry(a,tree,terminal,radial)
    a.out.mkdir(exist_ok=True,parents=True)

    d2=diameter_squared(cloud)
    scale=10**18
    t=F(isqrt(scale**2*d2.denominator//d2.numerator),scale)
    assert t*t*d2<=1<(t+F(1,scale))**2*d2
    f_max=t/(1-RHO)
    baseline = [sampled_group(n,arb.pi(),f_max*LAMBDA) for n in ORDERS]
    baseline.append(scale_sixth(cloud,f_max))
    unchanged=output_groups(a.out/"points-13668.json",baseline,f_max)
    assert unchanged["points_sha256"] == named("joint-contraction-audit.json")["point_sha256"]
    assert unchanged["total"]==13668
    variants={"max_area":([3009,579,348,249,193],10020),
              "smaller":([3007,579,348,248,193],10011)}
    receipts={}
    for name,(counts_by_arc,expected_total) in variants.items():
        pts,f,fits=tangent_variant(counts_by_arc,f_max,cloud)
        receipt=output_groups(a.out/f"points-{expected_total}.json",pts,f)
        assert receipt["total"]==expected_total and (f==f_max if name=="max_area" else f<f_max)
        receipt["tangent_cells"] = counts_by_arc
        receipt["direct_edge_enclosures"] = fits
        receipts[name]=receipt
    report=dict(target=str(TARGET), integer_area_floor=str(GOAL), shrink=str(LAMBDA), point_shrink=str(RHO),
                coverage=counts, intervals=geometry, baseline=unchanged, compressed=receipts,
                containment_count_limit=point_count_limit(),
                python=platform.python_version(), python_flint=flint.__version__, precision=ctx.prec,
                script_sha256=checksum(Path(__file__)))
    print(json.dumps(report,indent=2))

if __name__=="__main__":
    main()
