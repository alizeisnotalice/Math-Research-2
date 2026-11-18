"""Compare the two concrete continuous hosts, using the R215 true face segments."""

from fractions import Fraction as F
from pathlib import Path
import argparse, json, sys
from hashlib import sha256
from flint import arb, ctx


def materialize(ns, hs, hull, outer_rational):
    # Remove redundant nearly-parallel rows before enclosing adjacent intersections.
    # Dropping any row keeps the construction outward; positive angular coverage is
    # checked after the removals, so the validity of the upper bound is independent
    # of the choice of dropped lines.
    from collections import deque

    def intersection(i, j):
        rx, ry = ns[i]
        px, py = ns[j]
        det = rx * py - ry * px
        assert det > 0
        return hs[i] * py - ry * hs[j], rx * hs[j] - hs[i] * px, det

    def outside(k, i, j):
        x, y, d = intersection(i, j)
        return ns[k][0] * x + ns[k][1] * y > hs[k] * d

    active = deque()
    for k in range(len(ns)):
        while len(active) > 1 and outside(k, active[-2], active[-1]):
            active.pop()
        while len(active) > 1 and outside(k, active[0], active[1]):
            active.popleft()
        active.append(k)
    while len(active) > 2 and outside(active[0], active[-2], active[-1]):
        active.pop()
    while len(active) > 2 and outside(active[-1], active[0], active[1]):
        active.popleft()
    active = list(active)
    assert len(active) > 2
    corners = []
    for i, k in enumerate(active):
        x, y, d = intersection(active[i - 1], k)
        corners.append((F(x, d), F(y, d)))
    outer = outer_rational(corners, hull)
    return outer


def run(comp, repair, alt_in, pin_in):
    sys.path.insert(0, str(repair))
    sys.path.insert(0, str(comp))
    sys.path.insert(0, str(comp / "finite-b-core"))
    import certify_single_repair as R
    from certify_nonmonotone_ah_core import hull, double_area
    from deformed_reuleaux_support import body, environment
    from certify_ah_tapered_mirror_counterexample import Polygon, Support
    from certify_tapered_forced_support import Profile

    Q = R.Q
    ref = json.loads((repair / "single-repair.json").read_bytes())
    assert (
        ref["dependencies"]["verifier_sha256"]
        == sha256((repair / "certify_single_repair.py").read_bytes()).hexdigest()
    )
    receipts = [
        alt_in.with_name("continuous-153-area.json"),
        pin_in.with_name("single-H1-area.json"),
    ]
    records = [json.loads(z.read_bytes()) for z in receipts]
    assert all(
        z["config_sha256"] == sha256(path.read_bytes()).hexdigest()
        for path, z in zip([alt_in, pin_in], records)
    )
    old, core, container = R.refined_core(comp, 4096, hull)
    g = old["g"]
    ns = environment(str(comp), "1.5497")
    ctx.prec = 192
    folder = comp / "finite-b-core"
    orb = json.loads((folder / "host-poses.json").read_bytes())["poses"]

    def load(p):
        return json.loads(p.read_bytes())

    p = load(folder / "pentagons.json")
    c = load(folder / "curves.json")
    sample = hull(
        [(10 * x, 10 * y) for x, y in p["vertices_ccw_integer"] + c["combined_hull"]]
    )
    datasets = {s: load(comp / s) for s in set(row["source"] for row in orb)}
    mirror = Support(datasets["ah-tapered-mirror-counterexample-input.json"], arb.pi())
    mirrors = {"Support": mirror, "Polygon": Polygon(mirror)}
    old_h = []
    for rr in orb:
        if "w" in rr:
            vs, h = body(ns, F(rr["w"]), F(rr["e"]))
        elif "body_type" in rr:
            h = mirrors[rr["body_type"]].support_direction
        else:
            h = Profile(datasets[rr["source"]][rr["body_index"]]).support
        co, si = ns["trig"](R.exact(rr["turn"]))
        s = rr["reflection"]
        old_h.append((co, si, s, h))
    n, m, c3, b3, _, u, _, te, j = g.setup()
    half = arb(1) / 2

    def project(v, z):
        return g.plus(v, g.scale(half - g.dot(v, z), z))

    b1 = g.line_pair(n[0], n[1])
    f1 = g.line_pair(n[2], n[3])
    d1 = g.line_pair(n[4], n[5])
    e3 = g.line_pair(n[4], m[3])
    f3 = g.line_pair(n[3], m[2])
    G = g.circle_line(f3, n[5], 1)
    faces = [
        (c3, b1),
        (b1, project(e3, n[1])),
        (project(G, n[2]), f1),
        (f1, project(c3, n[3])),
        (e3, d1),
        (d1, G),
    ]
    for i, (a, b) in enumerate(faces):
        assert all((g.dot(n[i], x) - half).contains(0) for x in (a, b))
        assert all(h > g.dot(t, v) for _, t, h in old["otherEH"] for v in (a, b))
    oldbase = F(ref["base_area"][0])
    st = old["ns"]
    ta = g.plus(st["z"], st["D"])
    tp = old["q"](arb(1))
    chord = g.neg(g.rot(g.minus(tp, ta)))
    assert all(g.dot(chord, g.minus(v, ta)) < 0 for seg in faces for v in seg)
    pre = lambda v, co, si, s: (co * v[0] + s * si * v[1], -si * v[0] + s * co * v[1])

    def posed(w, e, k, s, tx, ty, pinned):
        V, h = body(ns, F(w), F(e))
        E = [g.minus(V[i], V[j]) for i, j in ns["orders"]]
        co, si = ns["trig"](R.exact(k))
        if pinned:
            g0 = half - h(pre(n[0], co, si, s))
            g2 = half - h(pre(n[2], co, si, s))
            tx, ty = g0, (2 * g2 + g0) / arb(3).sqrt()
        else:
            tx, ty = R.exact(tx), R.exact(ty)

        def hh(v):
            return h(pre(v, co, si, s)) + g.dot((tx, ty), v)

        def pp(v):
            z = pre(v, co, si, s)
            votes = [g.cross(a, z) for a in E]
            active = [
                i for i in range(10) if not votes[i] < 0 and not votes[(i + 1) % 10] > 0
            ]
            assert active
            out = []
            for i in active:
                j, r = ns["sector"][i]
                x, y = g.plus(V[j], g.scale(r, z))
                out.append((co * x - si * y + tx, s * (si * x + co * y) + ty))
            return out

        return hh, pp

    def boundbody(cfg):
        pin = "translation_rule" in cfg
        rows = [
            posed(*row, *t, pin)
            for row, t in zip(cfg["bodies"], cfg.get("translations", [["0", "0"]] * 3))
        ]
        public = core + sample
        seed = []
        # X samples are only candidate facet generators; every final support is evaluated exactly.
        for i in range(1800):
            nn = ns["trig"](R.exact(F(i, 5)))
            for _, p in rows:
                for x, y in p(nn):
                    seed.append((round(float(x.mid()) * Q), round(float(y.mid()) * Q)))
        P = hull(public + tuple(seed))
        normals = [(v - y, x - u) for (x, y), (u, v) in zip(P, P[1:] + P[:1])]
        aa, bb, cc = (
            R.support_scan(poly, normals) for poly in (core, sample, container)
        )
        rem = 12 if pin else 4
        rad = R.exact(F(rem, 8) * (F(".08464629") - F(".08464386")) ** 2)
        rr = R.exact(F(47, 10**10))
        assert 2 * (arb.pi() / (2 * 32768)).sin() ** 2 + R.exact(F(3, 10**11)) < rr
        hs = []
        count = 0

        def ceil(x):
            a, b = map(int, x.upper().man_exp())
            return a * 2**b if b >= 0 else -(-a // 2 ** (-b))

        for normal, ca, sa, pa in zip(normals, aa, bb, cc):
            x, y = normal
            ln = arb(x * x + y * y).sqrt()
            v = (arb(x) / ln, arb(y) / ln)
            B = ln * Q
            ob = min(pa, max(ca, sa + ceil(rr * B)))
            if ob > ca:
                # Only evaluate the 150 old bodies when the rounded samples could change the core facet.
                try:
                    fo = max(h(pre(v, co, si, s)).upper() for co, si, s, h in old_h)
                    ob = min(ob, max(ca, ceil(fo * B)))
                    count += 1
                except AssertionError:
                    pass
            fs = [h(v) for h, _ in rows]
            new = max(fs[0].upper(), (max(fs[1].upper(), fs[2].upper()) + rad).upper())
            hs.append(max(ob, ceil(new * B)))
        poly = materialize(normals, hs, hull, R.outer_rational)
        up = F(double_area(poly), 2 * Q**2)
        low = R.exact(oldbase)
        facesum = []
        for i, (a, b) in enumerate(faces):
            heights = [h(n[i]) for h, _ in rows]
            k = max(range(3), key=lambda j: float(heights[j].mid()))
            delta = heights[k] - half
            if not delta > 0:
                continue
            ps = rows[k][1](n[i])
            assert all(g.dot(v, n[j]) < half for v in ps for j in range(6) if j != i)
            inc = delta * g.norm(g.minus(b, a)) / 2
            facesum.append((i, str(inc)))
            low += inc
        ml, ex = map(int, low.lower().man_exp())
        lo = F(ml) * F(2) ** ex
        return dict(
            lower=str(lo),
            upper=str(up),
            faces=facesum,
            exact_old_directions=count,
            normal_count=len(normals),
        )

    alt = boundbody(load(alt_in))
    pin = boundbody(load(pin_in))
    lo = F(alt["lower"]) - F(pin["upper"])
    hi = F(alt["upper"]) - F(pin["lower"])
    assert F("3.36e-10") < lo < hi < F("1.074e-9")
    files = [
        repair / "certify_single_repair.py",
        repair / "single-repair.json",
        alt_in,
        pin_in,
        *receipts,
    ]
    return dict(
        alternating=alt,
        single=pin,
        difference=list(map(str, [lo, hi])),
        core_arc_subdivisions=4096,
        new_seed_directions=1800,
        input_sha256={j.name: sha256(j.read_bytes()).hexdigest() for j in files},
        script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for arg in ["comp", "repair", "alt", "pin", "out"]:
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    if not __debug__:
        raise RuntimeError("Assertions must be enabled.")
    d = run(*(x.resolve() for x in [a.comp, a.repair, a.alt, a.pin]))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(d, indent=2) + "\n")
    print("gap", *[float(F(q)) for q in d["difference"]])
