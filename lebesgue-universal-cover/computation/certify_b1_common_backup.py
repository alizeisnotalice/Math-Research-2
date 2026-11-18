"""One convex B1 candidate, its area budget, and a common exact backup pose near the R199 body."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import inspect
import json
from pathlib import Path

from flint import arb, ctx

import audit_optimized_eh_chord as eh_geom
from audit_eh_continuous_placement import polygon
from certify_ah_deep_single_line import green
from certify_ah_margin_envelope import replace_line
from certify_ah_three_line_cut import format_area
from certify_ah_two_hit_budgets import geometry
from certify_eh_upper_chain import front

PIVOT=F(".430")
CUT_ANGLE=F(".0187")
W=F(".429")
T=F("1e-8")
INPUT_RADIUS=F("1e-6")
FACET_SCALE=F("1e-9")
ANCHOR_OUTER_RADIUS=F("6e-10")
LOCAL_INPUT_RADIUS=F("1e-9")
SUBDIVISIONS=4096
W_HALF=F(".00001")
T_HALF=T
BASE_FILES=(
    "gibbs-1p5497-refined-18-128.json","eh-two-line-cover-certificate.json",
    "eh-upper-chain-parameters.json","flexible-006-eleven-eh.json",
    "ah-pivot-obstruction-certificate.json"
)


def rational_arb(x):
    f=F(x)
    return arb(f.numerator)/f.denominator


def bounded_area(archive,pivot,angle,bracket=".460"):
    text=replace_line(inspect.getsource(green),'F("0.07")',f'F("{bracket}")')
    impl=dict(green.__globals__)
    exec(compile(text,"b1_long_cap_bracket","exec"),impl)
    return impl["green"](archive,pivot,angle)


def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("archive",type=Path)
    args=ap.parse_args()
    archive=args.archive.resolve()
    raw={n:(archive/n).read_bytes() for n in BASE_FILES}
    refs={n:json.loads(value) for n,value in raw.items()}
    gibbs,old,params,baseline,obstruction=(refs[n] for n in BASE_FILES)
    assert F(gibbs["sigma_degrees"])==F("1.5497")
    assert F(old["summary"]["sigma_degrees"])==F("1.5497")
    assert baseline["R126_sha256"]==sha256(raw[BASE_FILES[1]]).hexdigest()
    assert baseline["parameter_sha256"]==sha256(raw[BASE_FILES[2]]).hexdigest()
    assert old["summary"]["AH_relocation_after_root"] and old["summary"]["delta_is_enclosing_placement_depth"]
    specimen=obstruction["explicit_witness"]
    assert F(specimen["reference_w"])==W and F(specimen["Reuleaux_triangle_Minkowski_weight"])==T
    assert F(specimen["pivot_w"])==PIVOT and F(specimen["AH_angle_radians"])==F(".006")
    ns,_,_=geometry(archive,[(F(".427"),F(".0048"))],far=F(".183"),negative=F(".210"),eps=F(".0051"))
    ctx.prec=160
    a=rational_arb
    pi=arb.pi()
    eh_geom.SIGMA=a(F("1.5497"))*eh_geom.PI/180
    eh=front(old,params)
    assert len(eh["active"])==11 and eh["active"]==baseline["active_facet_indices"]
    z,D,direct,normal=(ns[k] for k in ("z","D","direct","normal"))
    dot,add,sub,mul,perp,norm,cross=(ns[k] for k in ("dot","add","sub","mul","perp","norm","cross"))
    reflect=ns["reflection"]
    pivot=add(z,mul(a(PIVOT),D))
    n_ah=mul(-1,perp(add(mul(a(CUT_ANGLE).cos(),direct),mul(a(CUT_ANGLE).sin(),normal))))
    b_ah=dot(n_ah,pivot)
    # Both ears are disjoint. The complete modification of EH lies inside
    # the rectangle for the deepest enclosing EH chord.
    v,t=eh["v"],eh["t"]
    sep=min((b_ah-dot(n_ah,(x,y))).lower() for x in (v[0],t[0]) for y in (v[1],t[1]))
    assert sep>a(F(".1"))

    q=ns["centered_points"](W,W)
    edges=[sub(q[k],q[l]) for k,l in ns["orders"]]
    qa=ns["trig"](a(0))
    R=1/arb(3).sqrt()
    vert=((arb(0),R),(arb(1)/2,-R/2),(-arb(1)/2,-R/2))

    def tri_upper(u):
        zz=[dot(u,v) for v in vert]
        # Arc radius points are allowed when the (unit) outward direction
        # lies in the relevant normal sector. At ambiguous endpoints retain
        # that branch to preserve an upper bound.
        return max([z.upper() for z in zz]+[(1+z).upper() for z in zz if not z>arb("-0.5")])

    def support(u):
        return (1-a(T))*ns["supports"](q,edges,qa,reflect(u)).upper()+a(T)*tri_upper(u)

    paired=reflect(q[2])
    derivative=reflect(ns["point_jets"](a(W))[1][2])
    n_old=mul(-1/norm(derivative),perp(derivative))
    out=add(mul(1-a(T),paired),mul(a(T),vert[0]))
    violation=dot(n_old,sub(out,paired))
    assert violation>a(F("9.21748e-11"))
    n_hex=eh_geom.setup()[0]
    assert all(dot(n,out)<a(F(".5")) for n in n_hex)
    identity=add(mul(1-a(T),q[4]),mul(a(T),reflect(vert[0])))
    identity_violation=dot(n_ah,identity)-b_ah
    assert identity_violation>0

    # The deepest enclosing EH chord gives a smaller convex host than the
    # eleven-facet EH chain. Its inscribed chords, clipped by the new AH
    # line, form an inner polygon P of the actual curved candidate.
    inner,_=polygon(params[-1][0],SUBDIVISIONS)
    low_release,high_release=1744,1770
    n_release=[next(n for label,n,_ in inner if label==f"chord_Z_{j}")
               for j in (low_release,high_release)]
    assert cross(*n_release)>0 and dot(*n_release)>0
    assert cross(n_release[0],n_old)>0 and cross(n_old,n_release[1])>0
    for (_,n,_),(_,m,_) in zip(inner,inner[1:]+inner[:1],strict=True):
        assert cross(n,m)>0 and dot(n,m)>0
    inner.append(("AH",n_ah,b_ah))
    scaled={"name":None,"margin":arb(1)}
    local={"name":None,"margin":arb(1)}
    hex_facets={"B1-Q":1,"P-F1":2,"F1-H":3,"E3-D1":4,"D1-G":5,"C3-B1":0}
    added_edge_count=0
    for label,n,b in inner:
        assert b>a(F(".499")) and norm(n).contains(1)
        h=support(n)
        slack=(1+a(FACET_SCALE))*b-h
        assert slack>0,(label,str(slack))
        if slack<scaled["margin"]:
            scaled={"name":label,"margin":slack}
        added_edge_count+=1
        if label in hex_facets:
            assert b.contains(a(F(".5"))) and all(x.overlaps(y) for x,y in zip(n,n_hex[hex_facets[label]],strict=True))
            continue
        if label.startswith("chord_Z_") and low_release<int(label.removeprefix("chord_Z_"))<high_release:
            continue
        d=b-h
        assert d>a(LOCAL_INPUT_RADIUS),(label,str(d))
        if d<local["margin"]:
            local={"name":label,"margin":d}
    assert R*a(FACET_SCALE)<a(ANCHOR_OUTER_RADIUS)
    ah_gap=(b_ah-support(n_ah)).lower()
    assert ah_gap>a(INPUT_RADIUS) and ah_gap>a(F("1.08949e-6"))
    eh_rows=[(b-support(n)).lower() for n,b in zip(eh["normals"],eh["offsets"],strict=True)]
    assert min(eh_rows)>a(INPUT_RADIUS)

    # A fixed orientation is sufficient for the explicit moving-input box.
    _,prime=ns["point_jets"](ns["box"](W-W_HALF,W+W_HALF))
    derivative_upper=max(norm(v).upper() for v in prime)
    shape_shift=derivative_upper*a(W_HALF)+2*R*a(T_HALF)
    assert shape_shift<a(INPUT_RADIUS)

    lc,uc,shape=bounded_area(archive,PIVOT,CUT_ANGLE)
    lc0,uc0,_=bounded_area(archive,F(".427"),F(".006"),".203")
    gl,gu=map(F,gibbs["area_interval"])
    el,eu=map(F,old["summary"]["cut_gain_area_interval_exact"])
    gainl,gainu=map(F,baseline["EH_gain_over_R126_exact_interval"])
    assert eh["gain"].overlaps(a(gainl).union(a(gainu)))
    blo,bhi=gl-eu-gainu-uc,gu-el-gainl-lc
    oldlo,oldhi=map(F,baseline["area_exact_interval"])
    assert gl-eu-gainu-uc0==oldlo and gu-el-gainl-lc0==oldhi
    new_body_radius=INPUT_RADIUS+ANCHOR_OUTER_RADIUS
    perimeter_upper=F("3.464102")
    pi_upper=F("3.141593")
    assert a(perimeter_upper)>2*arb(3).sqrt() and a(pi_upper)>pi
    added=perimeter_upper*new_body_radius+pi_upper*new_body_radius**2
    small_added=perimeter_upper*ANCHOR_OUTER_RADIUS+pi_upper*ANCHOR_OUTER_RADIUS**2
    localized_radius=ANCHOR_OUTER_RADIUS+LOCAL_INPUT_RADIUS
    localized_added=perimeter_upper*localized_radius+pi_upper*localized_radius**2
    assert bhi+added<oldlo-F(".0000025")
    assert bhi+small_added<oldlo-F(".000006")
    assert bhi+localized_added<oldlo-F(".000006")
    assert blo<bhi

    print(json.dumps({
        "complete":True,"scope":"one fixed backup orientation without translation; specified inputs only",
        "source_hashes":{n:sha256(value).hexdigest() for n,value in raw.items()},
        "sigma_degrees":"1.5497","pivot":str(PIVOT),"AH_angle_radians":str(CUT_ANGLE),
        "w":str(W),"mix_weight":str(T),"pose":"reflect in the 60-degree axis; translation zero",
        "original_true_paired_tangent":{"normal":list(map(str,n_old)),"point":list(map(str,paired)),
            "released_point":list(map(str,out)),"violation_lower":str(violation.lower())},
        "deep_AH_identity_violation_lower":str(identity_violation.lower()),
        "rigid_anchor_AH_support_slack_lower":str(ah_gap),
        "rigid_anchor_EH_support_slack_lower":[str(v) for v in eh_rows],
        "anchor_AH_after_halo_lower":str((ah_gap-a(INPUT_RADIUS)).lower()),
        "anchor_EH_after_halo_lower":str((min(eh_rows)-a(INPUT_RADIUS)).lower()),
        "EH_rectangle_passes_AH_lower":str(sep),
        "full_inscribed_polygon_edge_count":added_edge_count,
        "max_required_dilation_factor_minus_one":str(FACET_SCALE),
        "all_dilated_support_rows_minimum":{"edge":scaled["name"],"margin":str(scaled["margin"])},
        "anchor_outward_Hausdorff_radius":str(ANCHOR_OUTER_RADIUS),
        "common_pose_input_Hausdorff_radius":str(INPUT_RADIUS),
        "wide_halo_global_outer_radius":str(new_body_radius),
        "local_release_cone":{"chord_indices":[low_release,high_release],
            "bounding_outward_normals":[list(map(str,n)) for n in n_release],
            "all_retained_nonhex_support_margin":{"edge":local["name"],"margin":str(local["margin"])},
            "common_pose_input_Hausdorff_radius":str(LOCAL_INPUT_RADIUS),
            "global_outer_radius":str(localized_radius),
            "added_area_upper":str(localized_added),"area_upper":str(bhi+localized_added),
            "area_upper_display":format_area(blo,bhi+localized_added)[1],
            "strict_saving_lower":str(oldlo-bhi-localized_added)},
        "explicit_input_family":{"w_range":[str(W-W_HALF),str(W+W_HALF)],
            "weight_range":[str(T-T_HALF),str(T+T_HALF)],"body_support_speed_w_upper":str(derivative_upper),
            "Hausdorff_from_anchor_upper":str(shape_shift.upper())},
        "base_area_interval":list(map(str,(blo,bhi))),"base_area_display":format_area(blo,bhi),
        "R199_area_interval":list(map(str,(oldlo,oldhi))),
        "AH_cap_interval":list(map(str,(lc,uc))),"AH_cut_green_integral":shape,
        "EH_total_cut_area_interval":list(map(str,(el+gainl,eu+gainu))),
        "local_convex_hull_area_addition_upper":str(added),
        "local_convex_hull_area_upper":str(bhi+added),
        "local_convex_hull_area_upper_display":format_area(blo,bhi+added)[1],
        "local_convex_hull_saved_area_lower":str(oldlo-bhi-added),
        "anchor_only_area_addition_upper":str(small_added),
        "anchor_only_area_upper":str(bhi+small_added),
        "anchor_only_area_upper_display":format_area(blo,bhi+small_added)[1],
        "anchor_only_saved_area_lower":str(oldlo-bhi-small_added)
    },indent=2))


if __name__=="__main__":
    main()
