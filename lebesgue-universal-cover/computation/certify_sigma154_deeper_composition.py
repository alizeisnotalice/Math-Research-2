"""Check a true 1.54-degree AH cap and compose all certified placements with EH."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from certify_three_sigma_area_budget import ah_clearance, exact, fixed_ah, pack, setup, subtract

SIGMA = "1.54"
D_CHART = F(".000555")
SINGLE = "sigma154-ah429-eh550-composition.json"
ELEVEN = "sigma154-ah429-eh-chain997-composition.json"
SINGLE_SHA = "c3244d638b51cec714d72a4012d807895d262e476c76ca03c52de84419bd3381"
ELEVEN_SHA = "d7a2c66a83d64f0afd64ea89ebee6e4ad6c4f1bde292eae286249d00ef8b12ee"
JET_SOURCES = (
    "certify_sigma_pinned_fan.py",
    "certify_sigma_jet_bracket.py",
    "certify_sigma_extended_graph.py",
)
SLOPE_SOURCES = (
    "certify_sigma_pinned_fan.py",
    "certify_sigma_derivative_brackets.py",
    "certify_sigma_extended_graph.py",
    "certify_sigma_derivative_classes.py",
)
RETREAT_SOURCE = "certify_sigma154_eh_row_retreat.py"


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_bytes())


def physical_chart(comp, cfg):
    from certify_sigma_extended_graph import geometry
    from certify_sigma_reachable_roots import reach

    pivot, turn = map(F, cfg["line"])
    far, negative, eps = map(F, cfg["physical_parameters"])
    assert F(cfg["sigma_degrees"]) == F(SIGMA)
    assert 0 < exact(turn) < arb.pi()/2 and 0 < far < 1 and negative > 0 and eps > 0
    ns, source, _ = geometry(comp, [(pivot, turn)], sigma=SIGMA, far=far, negative=negative, eps=eps)
    assert ns["cuts"][-1] == (pivot, turn)
    assert all(p < pivot and d < turn for p, d in ns["cuts"][:-1])
    assert F(ns["eh_depth"]) == D_CHART
    # Independent full-domain test, before any inclusion in a reference cell.
    A = ns["A"]
    left = ns["chart"](ns["add"](ns["z"], ns["mul"](A(ns["knots"][-1]), ns["D"])))[0]
    right = ns["coord_hull"](ns["psi"], left)
    gap = min(min((ell-left).lower(), (right-ell).lower())
              for pair in ns["parameter_boxes"] for ell in pair[:2])
    assert gap > 0

    # Rebuild the physical boxes to record the strict Q_A trigger margins.
    ys, ds, coordinates = [], [], []
    for _, _, _, a in ns["parameter_boxes"]:
        b = a + F(1, 5000)
        floor = ns["physical_floor"](a, b)
        height = ns["height_bound"](a)
        x = ns["add"](ns["add"](ns["z"], ns["mul"](ns["box"](a, b), ns["D"])),
                      ns["mul"](floor.union(height.upper()), ns["normal"]))
        y = x[1] - A(".518")
        u = ns["sub"](x, ns["gate_C"])
        d = ns["dot"](u, u) - 1
        assert y > 0 and d > 0
        ys.append(y.lower())
        ds.append(d.lower())
        coordinates.append(a)
    assert len(coordinates) == len(ns["physical_deficits"]) and coordinates
    bound = max(ns["deficit_by_bin"]).upper()
    assert bound < A(eps)
    ranges = {j: reach(ns, j) for j in range(len(ns["ref_boxes"]))}
    ranges = {j: [list(pair) for pair in rr] for j, rr in ranges.items() if rr}
    assert {j for j, d in enumerate(ns["deficit_by_bin"]) if not d.is_zero()} <= set(ranges)

    # Every true paired boundary point up to far has its negative projection
    # in this scan; the checked true cap intersection will be before far.
    total = (pivot + negative) * 5000
    assert total.denominator == 1 and ns["wp"] == pivot and ns["wf"] == far
    room = A(negative) - ns["physical_rights"][-1]
    assert room > 0
    row, = [(n, h) for label, n, h in ns["rows"] if label.startswith("U1-Te")]
    info = dict(
        geometry_sha256=sha256(source.encode()).hexdigest(),
        physical_rows=len(coordinates), physical_graph_endpoint=str(ns["wmax"]),
        reference_endpoint=str(ns["knots"][-1]), reference_gap=pack(gap),
        scan_w_endpoints=list(map(str, (-negative, pivot))),
        active_w_endpoints=list(map(str, (min(coordinates), max(coordinates) + F(1, 5000)))),
        negative_projection_room=pack(room), cap_box_height_upper=pack(ns["hmax"].upper()),
        first_hit_y_above_518=pack(min(ys)),
        first_hit_squared_distance_to_C_above_1=pack(min(ds)),
        deficit_upper=pack(bound), reference_ranges=ranges,
        reference_subcells=sum(b-a for rr in ranges.values() for a,b in rr),
    )
    return info, row


def area_geometry(comp, eh, cfg, chart_row):
    import audit_optimized_eh_chord as g

    single_path, eleven_path = eh / SINGLE, eh / ELEVEN
    assert sha(single_path) == SINGLE_SHA and sha(eleven_path) == ELEVEN_SHA
    single, eleven = load(single_path), load(eleven_path)
    assert eleven["single_sha256"] == SINGLE_SHA
    assert F(single["sigma"]) == F(eleven["sigma"]) == F(SIGMA)
    assert F(single["EH_chart_depth"]) == D_CHART
    assert arb(single["EH_root_upper"]) < 0 and arb(single["EH_placement_margin"]) > 0
    assert F(eleven["EH_right_critical_upper"]["exact"][1]) < 0 and arb(eleven["EH_pose_margin"]) > 0
    assert sorted(eleven["EH_11_active_facets"]) == list(range(11))
    assert F(single["chart_to_single_EH_nested_determinant"]["exact"][0]) > 0
    pivot, turn = map(F, cfg["line"])
    far = F(cfg["physical_parameters"][0])
    cap, crossing, slope = fixed_ah(comp, SIGMA, pivot, turn, far)
    assert 0 < crossing[0] < crossing[1] < far and slope > 0
    raw = subprocess.check_output([sys.executable, "-E", "-B", str(comp / "certify_gibbs_area_refined.py"),
                                   SIGMA, "16", "64"])
    base_file = json.loads(raw)
    assert base_file["passed"] and F(base_file["sigma_degrees"]) == F(SIGMA)
    base = tuple(map(F, base_file["area_interval"]))
    loose_base = tuple(map(F, single["base_area"]["exact"]))
    assert loose_base[0] <= base[0] < base[1] <= loose_base[1]
    geom = setup(SIGMA, g)
    U, C3, T = (geom[k] for k in ("U", "C3", "T"))
    radial = g.minus(U, C3)

    def end(d):
        c, s = g.unit(exact(d))
        return g.plus(C3, g.minus(g.scale(c, radial), g.scale(s, g.rot(radial))))

    def normal(v):
        direction = g.minus(T, v)
        return g.scale(-1 / g.norm(direction), g.rot(direction))

    chart_end = end(D_CHART)
    chart_n = normal(chart_end)
    row_n, row_h = chart_row
    assert all(a.overlaps(b) for a,b in zip(chart_n, row_n, strict=True))
    assert g.dot(chart_n, T).overlaps(row_h)
    results = []
    for label, d, rec in (
        ("single", F(single["EH_depth"]), single),
        ("eleven", F(eleven["EH_deepest"]), eleven),
    ):
        V = end(d)
        clearance = ah_clearance(geom, V, T, g, pivot, turn)
        gain = tuple(map(F, rec["EH_cut"]["exact"]))
        area = subtract(subtract(base, cap), gain)
        if label == "single":
            assert T[1] < V[1] < chart_end[1]
            assert g.cross(g.minus(T, chart_end), g.minus(V, chart_end)) < 0
            cost = F(0)
        else:
            assert T[1] < chart_end[1] < V[1]
            assert g.cross(g.minus(T, V), g.minus(chart_end, V)) < 0
            check = g.norm(g.minus(V, T)) * g.norm(g.minus(normal(V), chart_n))
            old = list(map(F, eleven["transfer_cost"]["exact"]))
            assert check.overlaps(exact(old[0]).union(exact(old[1])))
            cost = max(old[1], F(pack(check)["exact"][1]))
        results.append(dict(label=label, EH_deepest=str(d), AH_EH_clearance=pack(clearance),
                            EH_cut=pack(gain, 39), area=pack(area, 39), transfer_cost_upper=str(cost)))
    return dict(base=pack(base, 39), AH_cut=pack(cap, 39), AH_curve_crossing=pack(crossing, 35),
                AH_crossing_derivative_lower=str(slope), EH_first={p.name: sha(p) for p in (single_path,eleven_path)},
                gibbs_area_output_sha256=sha256(raw).hexdigest(), cases=results)


def verify_certificates(comp, cfg, chart, paths, costs, replay_work=None):
    ranges = chart["reference_ranges"]
    chosen = {}
    inputs = {}
    sources = {}
    for index, path in enumerate(paths):
        record = load(path)
        settings = record["parameters"]
        for key in ("line", "physical_parameters"):
            assert list(map(F, settings[key])) == list(map(F, cfg[key]))
        assert F(settings["sigma_degrees"]) == F(SIGMA)
        assert ("bins" in settings) != ("slope_bins" in settings)
        if "slope_bins" in settings:
            recipe = SLOPE_SOURCES
            bins = [(j, "slope", params) for j, params in settings["slope_bins"]]
            assert record["geometry_sha256"] == chart["geometry_sha256"]
        else:
            recipe = JET_SOURCES
            bins = settings["bins"]
        expected_sources = {name: sha(comp / name) for name in recipe}
        retreat = F(settings.get("EH_margin_retreat", 0))
        runner = comp / recipe[1]
        if retreat:
            assert retreat > 0
            runner = HERE / RETREAT_SOURCE
            assert record["geometry_sha256"] == chart["geometry_sha256"]
            tight = record["EH_tightening"]
            assert tight["label"].startswith("U1-Te")
            assert F(tight["retreat"]) == retreat and F(tight["depth"]) == D_CHART
            expected_sources[RETREAT_SOURCE] = sha(runner)
        else:
            assert "EH_tightening" not in record
        assert record["bits"] >= 128 and record["source_hashes"] == expected_sources
        sources.update(expected_sources)
        assert path.name not in inputs
        inputs[path.name] = sha(path)
        if replay_work is not None:
            ip = replay_work / f"jet-input-{index}.json"
            op = replay_work / f"jet-output-{index}.json"
            ip.write_text(json.dumps(settings) + "\n")
            subprocess.run([sys.executable, "-E", "-B", str(runner), str(comp), str(ip), str(op)],
                           check=True, stdout=sys.stderr)
            fresh = load(op)
            assert fresh["parameters"] == settings and fresh["rows"] == record["rows"]
            assert fresh["source_hashes"] == expected_sources
        for row, (j, mode, params) in zip(record["rows"], bins, strict=True):
            assert j == row["bin"]
            rr = ranges.get(j, [])
            assert rr == row["reachable"], (path.name, j, rr, row["reachable"])
            if not row["positive"] or not rr:
                continue
            assert mode in ("split", "direct", "slope")
            classes = row["classes"]
            nc = 1 if mode == "direct" else 2
            assert len(classes) == nc * len(rr)
            for k, subrange in enumerate(rr):
                tests = classes[k*nc:(k+1)*nc]
                assert all(t["reference_bin"] == j and t["subcells"] == subrange for t in tests)
                if mode == "split":
                    l, anchor, other, r, tau = map(F, params)
                    assert l < anchor and other < r and tau > 0
                    for test, (a,b,sign) in zip(tests, ((l,anchor,-1),(other,r,1)), strict=True):
                        assert list(map(F, test["bracket"])) == [a,b]
                        assert test["anchor_comparison_sign"] == sign
                        assert F(test["anchor_degrees"]) == anchor and F(test["anchor_threshold"]) == -tau
                elif mode == "slope":
                    cut, left_lo, right_lo, left_hi, right_hi = map(F, params)
                    assert left_lo < right_lo and left_hi < right_hi
                    for test, (a,b,sign) in zip(tests, ((left_lo,right_lo,-1),(left_hi,right_hi,1)), strict=True):
                        assert list(map(F, test["bracket"])) == [a,b]
                        condition = test["initial_slope_class"]
                        assert condition["comparison"] == sign and F(condition["cut"]) == cut
                else:
                    l, r = map(F, params[:2])
                    assert l < r and list(map(F, tests[0]["bracket"])) == [l,r]
                    assert tests[0]["pose"] == params[2:]
            root = min(arb(c["root_margin_lower"]).lower() for c in classes)
            pose = min(arb(c["placement_margin_lower"]).lower() for c in classes)
            assert root > 0 and pose > 0
            chosen.setdefault(j, []).append((pose, root, path.name, retreat))
    missing = sorted(set(ranges)-set(chosen))
    assert not missing, f"Missing verified AH bins {missing}"
    profiles = {}
    for label, cost in costs.items():
        options = {}
        for j, group in sorted(chosen.items()):
            def retained(record):
                pose, _, _, retreat = record
                return F(pack(pose)["exact"][0]) - max(cost-retreat, 0)

            best = max(group, key=retained)
            pose, root, name, retreat = best
            options[j] = dict(pose=pack(pose), sign=pack(root), certificate=name,
                              retreat=str(retreat), retained=str(retained(best)))
        worst, j = min((F(q["retained"]), j) for j,q in options.items())
        root, root_bin = min((F(q["sign"]["exact"][0]), n) for n,q in options.items())
        profiles[label] = dict(minimum_retained=str(worst), limiting_bin=j,
                               minimum_root=str(root), limiting_root_bin=root_bin, bins=options)
    return dict(
        certificate_sha256=inputs, source_sha256=sources, replayed=bool(replay_work),
        placement_profiles=profiles,
    )


def run(args):
    comp = args.comp.resolve()
    sys.path.insert(0, str(comp))
    cfg = load(args.config)
    ctx.prec = 160
    chart, row = physical_chart(comp, cfg)
    ctx.prec = 160
    measure = area_geometry(comp, args.eh.resolve(), cfg, row)
    replay_work = args.replay_inputs
    if replay_work is not None:
        replay_work.mkdir(parents=True, exist_ok=True)
        assert args.receipts
    decision = "true area and physical gate"
    if args.receipts:
        costs = {c["label"]: F(c["transfer_cost_upper"]) for c in measure["cases"]}
        frozen = verify_certificates(comp, cfg, chart, args.receipts, costs, replay_work)
        for case in measure["cases"]:
            selection = frozen["placement_profiles"][case["label"]]
            room = F(selection["minimum_retained"])
            case["retained_margin_lower"] = str(room)
            case["limiting_bin"] = selection["limiting_bin"]
            case["universal_cover"] = room > 0
        decision = "all first placements and AH double-hit bins verified"
    else:
        frozen = None
    return dict(parameters=cfg, decision=decision, chart=chart, geometry=measure, placement=frozen,
                config_sha256=sha(args.config), source_sha256=sha(Path(__file__)))


if __name__ == "__main__":
    assert __debug__
    parser = ArgumentParser(description=__doc__)
    for key in ("comp", "config"):
        parser.add_argument(key, type=Path)
    parser.add_argument("receipts", nargs="*", type=Path)
    parser.add_argument("--eh", type=Path, default=HERE)
    parser.add_argument("--replay-inputs", type=Path,
                        help="directory for parametrized independent AH proof replays")
    print(json.dumps(run(parser.parse_args()), separators=(",", ":")))
