"""Verify the geometric transfer of the R136 arithmetic to six rational groups."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import platform
from pathlib import Path

import flint
import numpy as np

from certify_finite_point_leaf import point_tables
from certify_nine import S, ar, append_point, arb, arc_body, ctx
from certify_tomar_six_high import HERE, Q, SourceTree, digest
from replay_tomar_083731_scored import kernel
from select_tomar_083731_scores import AREA, FINE, radius_rows, score_floor

TARGET = F(".83731")
GOAL = int(TARGET*AREA)+1
LAMBDA = F(".999999035")
RADIUS = F(".707")
INRADIUS = F(".42")
ORDERS = (0, 3, 5, 7, 9)


def lower_radius(n, pi):
    return ar(F(1,2)) if n == 0 else 1 - 1 / (2 * (pi/(2*n)).cos())


def sampled_group(n, pi):
    inradius = lower_radius(n, pi)
    # Half of the largest gap in centered unit-circle coordinates.
    def contained(k):
        half_gap = pi / (2*n*k) if n else pi / k
        sagitta = (1 - half_gap.cos()) * (1 if n else ar(F(1,2)))
        return sagitta + ar(F(12,S)) < ar(1-LAMBDA)*inradius

    k = next(k for k in range(1, 4096) if contained(k))
    points = []
    if n:
        arc_body(points, k, pi, n, arb(0), arb(0), arb(0))
    else:
        for j in range(k):
            sine, cosine = (2*pi*j/k).sin_cos()
            append_point(points, cosine/2, sine/2)
    points = sorted(set(points))
    assert len(points) == max(n,1)*k
    return {"n": n, "samples_per_arc": k, "denominator": S+10, "points": points}


def verified_table(path, tree, subdiv, commit):
    raw = path.read_bytes()
    if subdiv == 32:
        assert digest(raw) == commit["original_npz_sha256"]
    with np.load(path, allow_pickle=False) as tab:
        counts = tab["counts"]
        points = tab["points"]
        assert tab["grid"] == 1024 and tab["scale"] == Q
        assert str(tab["manifest_sha256"]) == digest((tree.path/"SHA256SUMS").read_bytes())
    assert counts.tolist() == [96, *(n*subdiv for n in ORDERS[1:])]
    assert points.shape == (1025, int(sum(counts)), 2) and points.dtype==np.int32
    stamp = sha256()
    for values in (counts, points):
        stamp.update(values.tobytes())
    assert stamp.hexdigest() == commit["canonical_arrays_sha256"]
    return counts, points, digest(raw)


def verify_core_table(tree, counts, points, pi):
    offsets = np.r_[0, np.cumsum(counts)]
    for b, n in enumerate(ORDERS):
        data = points[:, offsets[b]:offsets[b+1]].astype(np.int64)
        # Point-to-center integers fit safely in signed 64-bit arithmetic.
        assert np.max(np.abs(data)) < Q
        if not n:
            assert np.all(data == data[0])
            assert np.all((data[0]**2).sum(axis=1) < int(LAMBDA*Q/2)**2)
            continue
        anchor = F(0) if n == 3 else F(tree.root["hi"][3*b-4])
        center_radius = ar(LAMBDA)/(2*(pi/(2*n)).cos())
        radius_floor = int(LAMBDA*Q)
        if n == 3:
            assert np.all(data == data[0])
        for row in range(1 if n == 3 else 1025):
            for k in range(n):
                t = ar(anchor*row/1024) + 2*pi*k/n
                sin, cos = t.sin_cos()
                x,y = (Q*center_radius*value for value in (cos, sin))
                xx,yy = (round(float(z.mid())) for z in (x,y))
                assert abs(x-xx)<1 and abs(y-yy)<1 and max(abs(xx),abs(yy))<Q
                distance = np.abs(data[row] - np.array([xx,yy], dtype=np.int64))+1
                assert np.all(np.square(distance).sum(axis=1) < radius_floor**2)


def verify_boxes_and_caps(tree, pi):
    exterior = []
    rd = ar(LAMBDA/2)
    for i, n in enumerate(ORDERS[1:]):
        j = 3*i
        dist = min(F(tree.root["hi"][j]), F(tree.root["hi"][j+1]),
                   -F(tree.root["lo"][j]), -F(tree.root["lo"][j+1]))
        rn = ar(LAMBDA)*lower_radius(n, pi)
        lower = pi*(rd**2+rn**2)/2 + ar(dist)*(rd+rn)
        assert lower > ar(F(GOAL, AREA))
        exterior.append(str(lower))
    R = ar(RADIUS)
    radius_exterior = (pi-(rd/R).acos())*rd**2 + rd*(R**2-rd**2).sqrt()
    assert radius_exterior > arb(".838")
    assert 2*rd**2 > R**2
    # Compare both cap areas at the same point, moving the support wall inwards.
    assert all(lower_radius(n, pi) > ar(INRADIUS) for n in ORDERS)
    upper_factor = (1+LAMBDA)/4
    lower_factor = 4*INRADIUS*LAMBDA**2 - (1+LAMBDA)*(RADIUS-F(1,2))
    assert lower_factor > upper_factor
    return exterior, radius_exterior, upper_factor, lower_factor


def prior_proofs(a, tree):
    def read(name):
        raw = (HERE/name).read_bytes()
        return json.loads(raw), digest(raw)

    stage, stage_sha = read("tomar-083731-stage-ledger.json")
    assert stage["complete"] and stage["parts"] == 1024 and F(stage["target"]) == TARGET
    assert stage["integer_goal"] == GOAL
    seed, seed_sha = read("tomar-083731-seed-audit.json")
    assert stage["seed_receipt_sha"] == seed_sha and seed["target_integer"] == GOAL
    assert F(seed["target"]) == TARGET and F(stage["radial_target"]) == F(seed["radial_target"])
    assert seed["tree_manifest_sha256"] == digest((tree.path/"SHA256SUMS").read_bytes())
    seed_counts = seed["totals"]
    assert [seed_counts[k] for k in ("canonical","score_cover","radial_extra","remaining")] == [
        6710314, 5647793, 4662, 1057859
    ]
    assert seed_counts["canonical"] == sum(seed_counts[k] for k in ("score_cover", "radial_extra", "remaining"))
    assert seed["min_bucket"] == 2101 and score_floor(2101) == GOAL
    assert score_floor(2100) < GOAL
    radius_target = F(seed["radial_target"])
    assert radius_target == F(".83732")
    radial, radial_sha = read("five-body-radial-affine.json")
    assert seed["provenance"]["radial_certificate"] == radial_sha
    assert radius_rows(radial, radius_target) == seed["radii"]
    for row in seed["radii"]:
        assert F(row["integer_area_lower"], AREA) > radius_target
    for b, row in enumerate(seed["radii"]):
        for j in (3*b, 3*b+1):
            radius = F(row["integer_radius"], FINE)
            assert F(tree.root["lo"][j]) < -radius < radius < F(tree.root["hi"][j])

    # The four layers only use the published fan kernel, its scoring extension,
    # and the above certified radial integer floors. Rebuild its exact source.
    _, source_sha, _ = kernel(HERE, a.work, a.clang_cache)
    assert source_sha == stage["source_sha"]
    remain = seed_counts["remaining"]
    for layer, expected in zip(stage["audits"], [117228, 10146, 2011, 614], strict=True):
        assert layer["tested"] == remain == layer["accepted"]+layer["remaining"]
        assert layer["remaining"] == expected
        remain = expected

    leaves, leaves_sha = read("tomar-083731-retained.json")
    previous, previous_sha = read("finite-point-lift-083731.json")
    terminal, terminal_sha = read("finite-point-lift-083731-replay.json")
    assert previous["outer_sha256"] == terminal["outer_sha256"] == leaves_sha
    assert previous["outer_audit_sha256"] == terminal["outer_audit_sha256"] == stage_sha
    assert terminal["lift_sha256"] == previous_sha and terminal["complete"]
    assert F(terminal["target"]) == F(previous["target"]) == TARGET
    assert terminal["tree_manifest"] == seed["tree_manifest_sha256"]
    assert terminal["kernel_sha256"] == digest((HERE/"finite_point_leaf_guard.cpp").read_bytes())
    assert terminal["driver_sha256"] == digest((HERE/"replay_finite_point_global.py").read_bytes())
    keys = sorted((int(part), leaf) for part, ids in leaves.items() for leaf in ids)
    assert len(set(keys)) == len(keys) == remain == terminal["original_leaves"]
    assert [tuple(row[:2]) for row in previous["rows"]] == keys
    counts = terminal["successful_tree_totals"]
    assert counts == [sum(row[j] for row in previous["rows"]) for j in (5,6,7)]
    assert counts[0] == 2*(counts[1]+counts[2])-remain

    shape, shape_sha = read("finite-point-leaf-input.json")
    bank, bank_sha = read("finite-point-global-input.json")
    assert bank_sha == previous["prior_bank_sha256"] == terminal["bank_and_rows"]
    assert shape_sha == terminal["geometry_sha256"]
    full_shape = dict(shape, cells=bank["cells"])
    diameter_num, ptr, triples = point_tables(full_shape)
    assert diameter_num == terminal["diameter_squared_numerator"]
    assert len(ptr)-1 == terminal["cells"] == 384
    assert len(triples)//8 == terminal["balanced_triples"] == 2493
    assert len(shape["points"]) == terminal["points"] == 30
    sixth = dict(n="P30", denominator=shape["denominator"], points=shape["points"])
    return terminal, radial, seed, sixth, {
        "score_seed": seed_sha, "outer_stage": stage_sha, "leaves": leaves_sha, "terminal": terminal_sha,
        "sourced_outer_kernel": source_sha, "P30": shape_sha,
    }


def diameter_squared(group):
    pts = group["points"]
    maximum = max((x-u)**2+(y-v)**2 for j,(x,y) in enumerate(pts) for u,v in pts[:j])
    return F(maximum, group["denominator"]**2)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("tree", "table32", "table64", "work", "clang_cache"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    a.clang_cache.mkdir(parents=True, exist_ok=True)
    ctx.prec=160
    pi = arb.pi()
    tree = SourceTree(a.tree)
    terminal, radial_spec, seed, sixth, bindings = prior_proofs(a, tree)
    spec = json.loads((HERE/"tomar-six-high-input.json").read_bytes())
    old_grouped_raw = (HERE/"tomar-grouped-08371-ledger.json").read_bytes()
    outer = json.loads(old_grouped_raw)
    assert digest(old_grouped_raw) == seed["provenance"]["score_ledger"]
    assert not outer["historical_area_inputs"]
    old_index = (HERE/"tomar-grouped-08371-index.json").read_bytes()
    assert json.loads(old_index)["audit_sha256"] == digest(old_grouped_raw)
    assert digest(old_index) == seed["provenance"]["score_index"]
    radial_raw = (HERE/"five-body-radial-material.npz").read_bytes()
    assert digest(radial_raw) == radial_spec["table_sha"]
    assert digest((HERE/"certify_radial_affine_fans.py").read_bytes()) == radial_spec["source_sha256"]
    table_hashes = {}
    for subdiv in (32,64):
        ct, pts, stamp = verified_table(getattr(a, f"table{subdiv}"), tree, subdiv, spec["tables"][str(subdiv)])
        assert stamp == terminal[f"table{subdiv}"]
        verify_core_table(tree, ct, pts, pi)
        if subdiv == 32:
            assert stamp == outer["table_sha256"] == radial_spec["source_table_sha256"]
            with np.load(HERE/"five-body-radial-material.npz", allow_pickle=False) as radial:
                assert np.array_equal(radial["counts"], ct)
                assert np.array_equal(radial["points"], pts[0])
        table_hashes[str(subdiv)] = stamp
    exterior, radius_exterior, hi, lo = verify_boxes_and_caps(tree, pi)
    groups = [sampled_group(n, pi) for n in ORDERS]+[sixth]
    diameters = [diameter_squared(g) for g in groups]
    assert all(d<1 for d in diameters)
    raw = (json.dumps(groups, separators=(",", ":"))+"\n").encode()
    (a.work/"polygon-083731-points.json").write_bytes(raw)
    output = {
        "area_floor": str(F(GOAL,AREA)), "integer_floor": GOAL, "area_scale": AREA,
        "target": str(TARGET), "strict_margin": str(F(GOAL,AREA)-TARGET),
        "point_counts": [len(g["points"]) for g in groups],
        "total_points": sum(len(g["points"]) for g in groups),
        "core_scale": str(LAMBDA), "points_sha256": digest(raw),
        "group_diameter_squared": [str(r) for r in diameters],
        "tree_manifest_sha256": digest((tree.path/"SHA256SUMS").read_bytes()),
        "material_table_hashes": table_hashes,
        "proof_bindings": bindings,
        "all_material_coordinates_in_reduced_cores": True,
        "radial_table_is_row_zero": True,
        "root_exterior_area_lower": exterior,
        "radial_floors": [row["integer_area_lower"] for row in seed["radii"]],
        "area_lower_for_point_at_radius_.707": str(radius_exterior),
        "inward_wall_distance": str(INRADIUS*(1-LAMBDA)),
        "cap_square_comparison_factors": [str(hi),str(lo)],
        "python": platform.python_version(), "python_flint": flint.__version__,
    }
    print(json.dumps(output,indent=2))


if __name__=="__main__":
    main()
