"""Enclose the convex hull of the core and specified finite or interval families under relaxed strip constraints."""

import argparse, json, sys
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from flint import arb, arb_series, ctx


def track_dependencies(comp, spec):
    """Index every repository file the verifier opens, including dynamically read Python text."""
    import importlib.util

    project = comp.resolve().parent
    accessed = set()

    def track(event, args):
        if event != "open" or not isinstance(args[0], (str, bytes)):
            return
        name = args[0]
        if isinstance(name, bytes):
            name = name.decode()
        p = Path(name)
        if p.suffix == ".pyc":
            try:
                p = Path(importlib.util.source_from_cache(str(p)))
            except ValueError:
                pass
        p = p.resolve()
        if (
            p not in (Path(__file__).resolve(), spec.resolve())
            and p.is_relative_to(project)
            and p.is_file()
        ):
            accessed.add(p)

    sys.addaudithook(track)

    def manifest():
        import flint, platform

        paths = sorted(accessed)
        return dict(
            repository_files={
                str(p.relative_to(project)): sha256(p.read_bytes()).hexdigest()
                for p in paths
            },
            verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
            python=platform.python_version(),
            python_flint=flint.__version__,
        )

    return manifest


def run(comp, spec):
    comp = comp.resolve()
    sys.path.insert(0, str(comp))
    sys.path.insert(0, str(comp / "finite-b-core"))
    from certify_host_area import prove, hull, double_area, endpoints, TARGET
    from deformed_reuleaux_support import body, environment

    names = ["host-poses.json", "pentagons.json", "curves.json", "core-outer.json"]
    files = [comp / "finite-b-core" / f for f in names]
    baseline = prove(*files)
    old_poses = json.loads((comp / "finite-b-core" / "pose-audit.json").read_bytes())
    old_caps = json.loads((comp / "finite-b-core" / "obstruction.json").read_bytes())
    assert (
        old_poses["complete"]
        and old_poses["all_bodies_exactly_touch_all_six_hex_faces"]
        and old_poses["body_count"] == 150
    )
    assert (
        old_poses["poses_sha256"]
        == sha256((comp / "finite-b-core" / "host-poses.json").read_bytes()).hexdigest()
    )
    assert old_caps["complete"] and all(
        arb(v) <= arb(1) / 2 for v in old_caps["support_upper_150_plus_core"][2:4]
    )
    old_receipt = json.loads((comp / "finite-b-core" / "host-area.json").read_bytes())
    assert (
        baseline["host_area_upper"] == old_receipt["host_area_upper"]
        and baseline["input_sha256"] == old_receipt["input_sha256"]
    )
    _, ps, cv, core = [json.loads(f.read_bytes()) for f in files]
    V = core["vertex_grid"]
    n = ps["directions"]
    sc = V // ps["grid"]
    assert sc * ps["grid"] == V and n == cv["N"]
    points = list(map(tuple, core["outer_integer_vertices"])) + [
        (sc * x, sc * y) for x, y in ps["vertices_ccw_integer"] + cv["combined_hull"]
    ]
    ns = environment(str(comp), "1.5497")
    ctx.prec = 192
    A, dot, cross, rot = ns["A"], ns["dot"], ns["cross"], ns["perp"]
    config = json.loads(spec.read_bytes())
    single_body = []
    pinned = config.get("translation_rule") == "pin_hexagon_faces_0_and_2"
    if pinned:
        assert "translations" not in config
        # R215 proves that the old core and poses obey both original oblique caps.
        cap_directions = [ns["av"](ns["d"]["m"][i]) for i in (3, 5)]
    translations = config.get("translations", [["0", "0"] for _ in config["bodies"]])
    assert len(translations) == len(config["bodies"])
    for (w, e, k, s), (tx, ty) in zip(config["bodies"], translations, strict=True):
        tx, ty = map(lambda x: A(F(x)), (tx, ty))
        verts, support = body(ns, F(w), F(e))
        edges = [ns["sub"](verts[i], verts[j]) for i, j in ns["orders"]]
        co, si = ns["trig"](A(F(k)))

        def raw_support(nx, ny):
            return support((co * nx + s * si * ny, -si * nx + s * co * ny))

        if pinned:
            g0 = arb(1) / 2 - raw_support(*ns["trig"](A(0)))
            g2 = arb(1) / 2 - raw_support(*ns["trig"](A(120)))
            tx, ty = g0, (2 * g2 + g0) / arb(3).sqrt()
        hs = []
        for phi in (0, 60, 120, 180, 240, 300):
            nx, ny = ns["trig"](A(phi))
            h = raw_support(nx, ny) + tx * nx + ty * ny
            hs.append(str(h))
        rec = dict(
            w=w,
            e=e,
            turn=k,
            reflection=s,
            six_supports=hs,
            translation=list(map(str, (tx, ty))),
        )
        if pinned:
            # Exact equalities follow from the translation equation and width one.
            assert all(arb(hs[j]).contains(arb(1) / 2) for j in (0, 2, 3, 5))
            assert arb(hs[1]) > arb(1) / 2 and arb(hs[4]) < arb(1) / 2
            cap_heights = [
                raw_support(x, y) + tx * x + ty * y for x, y in cap_directions
            ]
            assert all(h < arb(1) / 2 for h in cap_heights)
            rec["cap_supports"] = list(map(str, cap_heights))
        single_body.append(rec)
        for j in range(n):
            nx, ny = ns["trig"](A(F(360 * j, n)))
            z = (co * nx + s * si * ny, -si * nx + s * co * ny)
            votes = [cross(ed, z) for ed in edges]
            sectors = [
                i for i in range(10) if not votes[i] < 0 and not votes[(i + 1) % 10] > 0
            ]
            assert sectors
            pair = None
            for l in sectors:
                i, r = ns["sector"][l]
                raw = ns["add"](verts[i], ns["mul"](r, z))
                p = (
                    co * raw[0] - si * raw[1] + tx,
                    s * (si * raw[0] + co * raw[1]) + ty,
                )
                if pair is None:
                    pair = tuple(round(float((a * V).mid())) for a in p)
                assert all(abs(v * V - q) < 1 for v, q in zip(p, pair))
            points.append(pair)
    pts = hull(points)
    assert pts
    extra_radius = arb(0)
    interval_audit = None
    if "family_interval" in config:
        inds = config["family_interval"]
        left = config["bodies"][inds["first"]]
        right = config["bodies"][inds["last"]]
        assert left[0] == right[0] and left[2:] == right[2:]
        w, l, r = F(left[0]), F(left[1]), F(right[1])
        assert l < r
        ebox = ns["box"](l, r)
        body(ns, w, ebox)
        add, mul = ns["add"], ns["mul"]
        ell, t0 = ns["chart"](add(ns["z"], mul(A(w), ns["D"])))
        p = ns["psi"](ell)
        t = t0 - arb_series([ebox, 1])
        # Two anchors are constant and the middle contact is affine.
        # The reflection of X(p,-t) is an isometry.
        changing = [ns["X"](ell, t), ns["X"](p, -t)]
        second = [2 * (x[2] ** 2 + y[2] ** 2).sqrt() for x, y in changing]
        assert max(second) < A(4)
        # Interpolating the exact pinned translations adds at most 8 to the
        # acceleration. The inverse two-normal matrix has operator norm sqrt(2).
        curvature_bound = 12 if pinned else 4
        extra_radius = A(curvature_bound) * A(r - l) ** 2 / 8
        if pinned:
            relevant = [single_body[inds[i]] for i in ("first", "last")]
            assert all(
                arb(q["six_supports"][1]) - extra_radius > arb(1) / 2 for q in relevant
            )
            assert all(
                arb(q["six_supports"][1]) + extra_radius < arb(".500000301")
                for q in relevant
            )
            assert all(
                arb(q["six_supports"][4]) + extra_radius < arb(1) / 2 for q in relevant
            )
            assert all(
                arb(h) + extra_radius < arb(1) / 2
                for q in relevant
                for h in q["cap_supports"]
            )
        interval_audit = dict(
            w=str(w),
            ends=list(map(str, [l, r])),
            anchor_acceleration_bounds=list(map(str, second)),
            interpolation_radius=str(extra_radius),
            radius_formula=f"{curvature_bound}*(r-l)^2/8",
        )

    def upper(poly):
        area = F(double_area(poly), 2 * V**2)
        perimeter = (
            sum(
                (
                    arb((x - u) ** 2 + (y - v) ** 2).sqrt()
                    for (x, y), (u, v) in zip(poly, poly[1:] + poly[:1])
                ),
                arb(0),
            )
            / V
        )
        r = 2 * (arb.pi() / (2 * n)).sin() ** 2 + A(F(3, 10**11)) + extra_radius
        return area + endpoints(r * perimeter + arb.pi() * r * r)[1]

    unrelaxed = upper(pts)
    print("area", float(unrelaxed), flush=True)
    assert unrelaxed < TARGET
    width_lower = []
    for i in range(3):
        values = [
            max(
                [arb(1) / 2]
                + [arb(row["six_supports"][j]).lower() for row in single_body]
            )
            for j in (i, i + 3)
        ]
        width_lower.append(str((values[0] + values[1]).lower()))
    return dict(
        N=n,
        interval_audit=interval_audit,
        width_lower=width_lower,
        six_body_supports=single_body,
        family=config["bodies"],
        translations=None if pinned else translations,
        unrelaxed_upper=str(unrelaxed),
        translation_rule=config.get(
            "translation_rule", "affine endpoint interpolation"
        ),
        B=str(TARGET),
        previous_outer_inputs=baseline["input_sha256"],
        original_150_upper=baseline["host_area_upper"],
        corrected_vertex_grid=V,
        config_sha256=sha256(spec.read_bytes()).hexdigest(),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("comp", type=Path)
    ap.add_argument("spec", type=Path)
    ap.add_argument("out", type=Path)
    a = ap.parse_args()
    if not __debug__:
        raise RuntimeError("Assertions must be enabled.")
    manifest = track_dependencies(a.comp, a.spec)
    r = run(a.comp, a.spec)
    r["dependencies"] = manifest()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(r, indent=2) + "\n")
    print(float(F(r["unrelaxed_upper"])))
