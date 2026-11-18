"""Necessary area after posing one further exact input in a convex superset of the fixed core."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys
from flint import arb, ctx


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
    from certify_host_area import prove

    baseline_paths = [
        comp / "finite-b-core" / name
        for name in (
            "host-poses.json",
            "pentagons.json",
            "curves.json",
            "core-outer.json",
        )
    ]
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
    rebuilt = prove(*baseline_paths)
    assert (
        rebuilt["host_area_upper"] == old_receipt["host_area_upper"]
        and rebuilt["input_sha256"] == old_receipt["input_sha256"]
    )
    from certify_ah_tapered_minkowski_closure import geometry
    from deformed_reuleaux_support import body, environment
    from certify_nonmonotone_ah_core import hull, double_area, BASE_LO, LENGTH2_LO

    old, _ = geometry(comp)
    ns = environment(str(comp), "1.5497")
    ctx.prec = 160
    A, g = ns["A"], old["g"]
    ta = g.plus(ns["z"], ns["D"])
    tp = old["q"](arb(1))
    ch = g.minus(tp, ta)
    lensq = g.dot(ch, ch)
    core = old["oldarea"] + old["oldcap"] - old["green"](F(0), F(1), ta)
    assert core > A(BASE_LO) and lensq > A(LENGTH2_LO)
    scale = 10**11
    err = arb(0).union(arb(2) / scale).union(-arb(2) / scale)

    def inner(p):
        d = g.minus(p, ta)
        x = g.dot(d, ch) / lensq
        y = -g.cross(ch, d) / lensq
        if not y > arb(10000) / scale:
            return None
        iy = int(float(y.mid()) * scale) - 800
        t = arb(iy) / (scale * y)
        assert 0 < t < 1
        left, right = t * x, t * x + 1 - t
        ix = round(float(((left + right) * scale / 2).mid()))
        if not scale * left < ix < scale * right:
            return None
        return ix, iy

    initial = {(0, 0), (scale, 0)}
    dig = {}
    for name, key, count in [
        ("pentagons.json", "vertices_ccw_integer", 117),
        ("curves.json", "combined_hull", 33),
    ]:
        path = comp / "finite-b-core" / name
        dat = json.loads(path.read_bytes())
        dig[name] = sha256(path.read_bytes()).hexdigest()
        if count == 117:
            assert (
                dat["body_count"] == 117
                and dat["grid"] == scale
                and F(dat["coordinate_rounding_error_strict"]) <= F(2, scale)
            )
        else:
            assert (
                dat["complete"]
                and len(dat["bodies"]) == 33
                and dat["scale"] == scale
                and F(dat["coordinate_error_asserted"]) <= F(2, scale)
            )
        for px, py in dat[key]:
            q = inner((arb(px) / scale + err, arb(py) / scale + err))
            if q is not None:
                initial.add(q)
    base = hull(initial)
    from certify_nonideal_pentagon_pose_library import verify

    def bound(poly):
        return BASE_LO + LENGTH2_LO * double_area(poly) / (2 * scale**2)

    n, m, c3, b3, *_ = g.setup()

    def projected(p, normal):
        return g.plus(p, g.scale(arb(1) / 2 - g.dot(p, normal), normal))

    e3 = g.line_pair(n[4], m[3])
    f3 = g.line_pair(n[3], m[2])
    segments = [(projected(b3, m[3]), e3), (projected(f3, m[5]), c3)]
    chordnormal = g.neg(g.rot(ch))
    for p, q in segments:
        assert g.norm(g.minus(p, q)) > A(F(".2535"))
        assert all(h > g.dot(t, v) for _, t, h in old["otherEH"] for v in (p, q))
        assert all(g.dot(chordnormal, g.minus(v, ta)) < 0 for v in (p, q))
    base_area = bound(base)
    B = F(".844082046521218581802948589770767704354")
    config = json.loads(spec.read_bytes())
    rows = []
    for w, e, turn, sgn in config["bodies"]:
        w, e = F(w), F(e)
        proof = verify(ns, w, e, F(0), F(0))
        caps = [arb(j[3]).lower() for j in proof["cap_eliminations"]]
        # These triangles have base on C150's true cap face, so outside C150.
        # Use only its fixed base area, even when later forced copies cross that cap.
        other = A(base_area) + A(F(".2535")) * min(caps) / 2
        assert other > A(B)
        vertices, support = body(ns, w, e)
        edges = [g.minus(vertices[i], vertices[j]) for i, j in ns["orders"]]
        options = {}
        for k, s in [(0, 1), (240, -1)]:
            pts = set(base)
            for i in range(1201):
                n = ns["trig"](A(75 + F(i, 12)))
                pre = old["preimage"](n, k, s)
                sectors = [
                    j
                    for j in range(10)
                    if g.cross(edges[j], pre) > 0
                    and g.cross(pre, edges[(j + 1) % 10]) > 0
                ]
                if not sectors:
                    continue
                idx, rad = ns["sector"][sectors[0]]
                h = g.dot(vertices[idx], pre) + rad
                deriv = g.dot(vertices[idx], g.rot(pre))
                v = g.plus(g.scale(h, n), g.scale(s * deriv, g.rot(n)))
                q = inner(v)
                if q is not None:
                    pts.add(q)
            options[(k, s)] = hull(pts)
        keep = (turn, sgn)
        reject = next(p for p in options if p != keep)
        assert bound(options[reject]) > B, (
            str(w),
            str(e),
            float(bound(options[reject])),
        )
        base = options[keep]
        r = dict(
            w=str(w),
            e=str(e),
            forced=list(keep),
            rejected=list(reject),
            lower_rejected=str(bound(options[reject])),
            lower_forced=str(bound(base)),
            min_other_ten_cap_gap=str(min(caps)),
            other_ten_area_lower=str(other.lower()),
            root_certificate=proof,
        )
        rows.append(r)
        print(w, e, keep, float(bound(options[reject])), float(bound(base)), flush=True)
    return dict(
        B=str(B),
        base_lower=str(base_area),
        source_digests=dig,
        config_sha256=sha256(spec.read_bytes()).hexdigest(),
        stages=rows,
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
    obj = run(a.comp, a.spec)
    obj["dependencies"] = manifest()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(obj, indent=2) + "\n")
