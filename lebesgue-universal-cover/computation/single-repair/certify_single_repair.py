"""Certify the cost of adding just the R215 obstruction, at both possible poses."""

from fractions import Fraction as F
import argparse, json, sys
from hashlib import sha256
from math import isqrt
from pathlib import Path
from flint import arb, ctx

Q = 10**12
N = 32768


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


def exact(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def thick_grid(p):
    # Every point in the Arb rectangle is enclosed by the four returned lattice vertices.
    a, b = p
    X = round(float((a * Q).mid()))
    Y = round(float((b * Q).mid()))
    dx = max(1, int(float(abs(a * Q - X).upper())) + 2)
    dy = max(1, int(float(abs(b * Q - Y).upper())) + 2)
    assert abs(a * Q - X) < dx and abs(b * Q - Y) < dy
    return [(X + i * dx, Y + j * dy) for i in (-1, 1) for j in (-1, 1)]


def clip(poly, a, b, h):
    out = []
    for i, (x, y) in enumerate(poly):
        u, v = poly[i - 1]
        A = a * u + b * v - h
        B = a * x + b * y - h
        if (A > 0) != (B > 0):
            t = F(A, A - B)
            out.append((u + t * (x - u), v + t * (y - v)))
        if B <= 0:
            out.append((x, y))
    assert out
    return out


def outer_rational(poly, hull):
    corners = []
    for x, y in poly:
        a = int(x // 1)
        b = int(y // 1)
        corners.extend([(a, b), (a + 1, b), (a, b + 1), (a + 1, b + 1)])
    return hull(corners)


def row(u, h):
    # Coordinates used by clip are scaled by Q; all true points lie in [-1,1]^2.
    a, b = (round(float((v * Q).mid())) for v in u)
    hi = h + sum(abs(x - arb(n) / Q) for x, n in zip(u, (a, b)))
    c = int(float((hi * Q).upper())) + 3
    while not arb(c) / Q > hi:
        c += 1
    return a, b, c * Q


def refined_core(comp, subdiv, hull):
    from certify_ah_tapered_minkowski_closure import geometry

    old, _ = geometry(comp)
    g = old["g"]
    ns = old["ns"]
    n, m, c3, b3, _, u, _, te, j = g.setup()
    pre = arb(1) / 2
    lip = g.line_pair

    def project(v, z):
        return g.plus(v, g.scale(pre - g.dot(v, z), z))

    b1 = lip(n[0], n[1])
    f1 = lip(n[2], n[3])
    d1 = lip(n[4], n[5])
    e3 = lip(n[4], m[3])
    f3 = lip(n[3], m[2])
    G = g.circle_line(f3, n[5], 1)
    q = project(e3, n[1])
    p = project(G, n[2])
    I = project(b3, m[3])
    H = project(c3, n[3])
    K = project(f3, m[5])
    c = (-e3[0], e3[1])
    Y = (arb(0), c[1] + (1 - c[0] ** 2).sqrt())
    ta = g.plus(ns["z"], ns["D"])
    tp = old["q"](arb(1))
    points = [b1, q, Y, ta, tp, p, f1, H, j, I, e3, d1, G, K, c3, b1]
    centres = {1: e3, 2: c, 4: G, 7: c3, 8: b3, 12: f3}
    outer = []
    for i in range(len(points) - 1):
        a, b = points[i : i + 2]
        outer.extend(thick_grid(a))
        if i not in centres:
            continue
        O = centres[i]
        r = g.minus(a, O)
        s = g.minus(b, O)
        assert g.cross(r, s) > 0 and g.dot(r, s) > 0
        assert g.dot(r, r).contains(1) and g.dot(s, s).contains(1)
        prev = r
        for k in range(1, subdiv + 1):
            t = arb(k) / subdiv
            v = g.plus(g.scale(1 - t, r), g.scale(t, s))
            z = g.scale(1 / g.norm(v), v)
            # Exact intersection of two adjacent unit-circle tangent lines.
            contact = g.plus(O, g.scale(1 / (1 + g.dot(prev, z)), g.plus(prev, z)))
            outer.extend(thick_grid(contact))
            prev = z
    poly = hull(outer)
    for _, v, h in old["otherEH"]:
        poly = clip(poly, *row(v, h))
    core = outer_rational(poly, hull)
    container = [(-Q, -Q), (Q, -Q), (Q, Q), (-Q, Q)]
    for v in [*n, m[3], m[5]]:
        container = clip(container, *row(v, pre))
    return old, core, outer_rational(container, hull)


def support_scan(vertices, normals):
    def q(v, u):
        return v[0] * u[0] + v[1] * u[1]

    L = len(vertices)
    k = max(range(L), key=lambda j: q(vertices[j], normals[0]))
    ans = []
    for u in normals:
        while q(vertices[(k + 1) % L], u) > q(vertices[k], u):
            k = (k + 1) % L
        # Ties are exact and walking beyond the last equal point cannot skip an extremum.
        while q(vertices[(k + 1) % L], u) == q(vertices[k], u):
            k = (k + 1) % L
        assert q(vertices[(k - 1) % L], u) <= q(vertices[k], u) and q(
            vertices[(k + 1) % L], u
        ) <= q(vertices[k], u)
        ans.append(q(vertices[k], u))
    return ans


def envelope(core, raw, container, hull, double_area):
    sample = hull(raw)
    P = hull(core + sample)
    ns = [(v - y, x - u) for (x, y), (u, v) in zip(P, P[1:] + P[:1])]
    a, b, c = (support_scan(poly, ns) for poly in (core, sample, container))
    # The width-one boundary has curvature radius between 0 and 1; its nearest angular sample
    # loses at most 1-cos(pi/N). All old sample point errors are below 3e-11 in norm.
    rad = 2 * (arb.pi() / (2 * N)).sin() ** 2 + exact(F(3, 10**11))
    dscale = 10**16
    rr = int(float((rad * dscale).upper())) + 2
    assert arb(rr) / dscale > rad
    hs = []
    for n, ca, sa, pa in zip(ns, a, b, c):
        norm = isqrt(n[0] ** 2 + n[1] ** 2) + 1
        shift = (rr * Q * norm + dscale - 1) // dscale
        hs.append(min(pa, max(ca, sa + shift)))
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
    return F(double_area(outer), 2 * Q**2), len(outer)


def sample_new(
    ns, w, e, pose, point_parameters=False, inner=None, common=None, stride=1
):
    co, si = ns["trig"](exact(pose[0]))
    s = pose[1]
    from deformed_reuleaux_support import body

    vs, h = body(ns, w, e)
    edges = [ns["sub"](vs[i], vs[j]) for i, j in ns["orders"]]
    # The normal fan pins the marked contact anchors. Their projections are identically half.
    for j, idx in enumerate((1, 3, 0)):
        n = ns["mul"](-1, ns["av"](ns["d"]["n"][j]))
        slots = [
            k
            for k in range(10)
            if ns["cross"](edges[k], n) > 0 and ns["cross"](n, edges[(k + 1) % 10]) > 0
        ]
        assert len(slots) == 1 and tuple(ns["sector"][slots[0]]) == (idx, 0)
        assert ns["dot"](vs[idx], n).contains(arb(1) / 2)
    for j in (3, 5):
        x, y = ns["av"](ns["d"]["m"][j])
        pre = (co * x + s * si * y, -si * x + s * co * y)
        assert h(pre) < arb(1) / 2
    out = []
    for j in range(0, N, stride):
        x, y = ns["trig"](exact(F(360 * j, N)))
        d = (co * x + s * si * y, -si * x + s * co * y)
        cv = [ns["cross"](v, d) for v in edges]
        active = [i for i in range(10) if not cv[i] < 0 and not cv[(i + 1) % 10] > 0]
        assert active
        union = None
        for i in active:
            k, r = ns["sector"][i]
            a, b = ns["add"](vs[k], ns["mul"](r, d))
            v = (co * a - si * b, s * (si * a + co * b))
            union = v if union is None else tuple(a.union(b) for a, b in zip(union, v))
            corners = thick_grid(v)
            if point_parameters:
                assert all(abs(v[j] * Q - c[j]) < 6 for c in corners for j in (0, 1))
            out.extend(corners)
        if inner is not None:
            witness = inner(union)
            if witness is not None:
                common.add(witness)
    return out


def run(comp, subdiv):
    comp = comp.resolve()
    sys.path.insert(0, str(comp))
    sys.path.insert(0, str(comp / "finite-b-core"))
    from certify_nonmonotone_ah_core import hull, double_area, BASE_LO, LENGTH2_LO
    from deformed_reuleaux_support import environment
    from certify_nonideal_pentagon_pose_library import verify
    from certify_host_area import prove

    names = ["host-poses.json", "pentagons.json", "curves.json", "core-outer.json"]
    paths = [comp / "finite-b-core" / name for name in names]
    baseline = prove(*paths)
    assert baseline["N"] == N
    certs = [json.loads(p.read_bytes()) for p in paths]
    capcert = json.loads((comp / "finite-b-core" / "obstruction.json").read_bytes())
    assert (
        capcert["complete"]
        and len(capcert["all_angle_certificate"]["cap_eliminations"]) == 10
    )
    ctx.prec = 192
    old, core, container = refined_core(comp, subdiv, hull)
    ns = environment(str(comp), "1.5497")
    ctx.prec = 192
    z = ns["A"]
    w, e, rw, re = map(F, (".7982907", ".032688", "1e-10", "1e-11"))
    rootcert = verify(ns, w, e, rw, re)
    raw = [
        (10 * x, 10 * y)
        for x, y in certs[1]["vertices_ccw_integer"] + certs[2]["combined_hull"]
    ]
    orig = hull(raw)
    up, _ = envelope(core, orig, container, hull, double_area)
    # Inscribed triangles supply independent area lower endpoints (inward treatment of rounded points).
    g = old["g"]
    st = old["ns"]
    ta = g.plus(st["z"], st["D"])
    tp = old["q"](arb(1))
    ch = g.minus(tp, ta)
    lensq = g.dot(ch, ch)
    corearea = old["oldarea"] + old["oldcap"] - old["green"](F(0), F(1), ta)
    assert corearea > z(BASE_LO) and lensq > z(LENGTH2_LO)
    S = 10**11
    err = arb(0).union(exact(F(2, S))).union(-exact(F(2, S)))

    def inner(p):
        d = g.minus(p, ta)
        x = g.dot(d, ch) / lensq
        y = -g.cross(ch, d) / lensq
        if not y > arb(10000) / S:
            return None
        for shrink in (2, 4, 8, 16, 32, 64, 128, 256, 512, 800, 1600):
            iy = int(float(y.mid()) * S) - shrink
            t = arb(iy) / (S * y)
            if not 0 < t < 1:
                continue
            l, r = t * x, t * x + 1 - t
            ix = round(float(((l + r) * S / 2).mid()))
            if S * l < ix < S * r:
                return ix, iy
        return None

    inside = {(0, 0), (S, 0)}
    for x, y in orig:
        m = inner((arb(x) / Q + err, arb(y) / Q + err))
        if m is not None:
            inside.add(m)

    def lower(poly):
        return BASE_LO + LENGTH2_LO * double_area(hull(poly)) / (2 * S**2)

    tol = F(1, 10**30)
    lo = lower(inside) - tol
    up += tol
    result = {}
    for k, s in [(0, 1), (240, -1)]:
        pts = sample_new(ns, w, e, (k, s), True)
        shot = hull(pts)
        end, vcount = envelope(core, hull(orig + shot), container, hull, double_area)
        # Each rounded corner of a point rectangle lies within 6/Q of the exact
        # support point (the present inputs are point parameters).
        newerr = arb(0).union(arb(6) / Q).union(-arb(6) / Q)
        add = set(inside)
        for x, y in shot:
            m = inner((arb(x) / Q + newerr, arb(y) / Q + newerr))
            if m is not None:
                add.add(m)
        start = lower(add) - tol
        end += tol
        result[f"{k},{s}"] = dict(
            area=list(map(str, (start, end))),
            increase=list(map(str, (start - up, end - lo))),
            outer_vertices=vcount,
        )
    box = (ns["box"](w - rw, w + rw), ns["box"](e - re, e + re))
    identity_common = set(inside)
    reflected_common = set(inside)
    sample_new(ns, *box, (0, 1), inner=inner, common=identity_common, stride=8)
    ptsbox = hull(sample_new(ns, *box, (240, -1), inner=inner, common=reflected_common))
    all_identity_lower = lower(identity_common) - tol
    all_reflected_lower = lower(reflected_common) - tol
    boxed, _ = envelope(core, hull(orig + ptsbox), container, hull, double_area)
    boxed += tol
    B = F(baseline["B"])
    previous_upper = F(baseline["host_area_upper"])
    previous_lower = F(
        json.loads((comp / "nonmonotone-ah-core.json").read_bytes())["guaranteed_lower"]
    )
    idlo, idhi = map(F, result["0,1"]["increase"])
    reflo, refhi = map(F, result["240,-1"]["increase"])
    assert refhi < B - previous_upper and idlo > B - previous_lower and boxed < B
    assert all_identity_lower > B and all_reflected_lower > lo
    assert refhi < idlo and idhi - idlo < F("3e-9") and refhi - reflo < F("3e-9")
    return dict(
        subdivisions=subdiv,
        directions=N,
        strict_endpoints=True,
        parameters={"w": str(w), "e": str(e), "w_radius": str(rw), "e_radius": str(re)},
        core_upper=str(F(double_area(core), 2 * Q**2)),
        base_area=list(map(str, (lo, up))),
        poses=result,
        each_box_identity_area_lower=str(all_identity_lower),
        each_box_reflected_area_lower=str(all_reflected_lower),
        box_quantifiers="One simultaneous host including every member. Identity individually exceeds B on the whole closed box.",
        box_reflected_area=[result["240,-1"]["area"][0], str(boxed)],
        box_reflected_increase=[str(reflo), str(boxed - lo)],
        B=str(B),
        box_reflected_slack_lower=str(B - boxed),
        R215_budget_conservative=str(B - previous_upper),
        old_family_budget=str(B - previous_lower),
        root_cert=rootcert,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("comp", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--subdivisions", type=int, default=2048)
    a = p.parse_args()
    if not __debug__:
        raise RuntimeError("Assertions must be enabled.")
    manifest = track_dependencies(a.comp, Path(__file__))
    r = run(a.comp, a.subdivisions)
    r["dependencies"] = manifest()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(r, indent=2) + "\n")
    print("base", list(map(lambda x: float(F(x)), r["base_area"])))
    for k, v in r["poses"].items():
        print(k, {l: [float(F(x)) for x in v[l]] for l in ("area", "increase")})
    print("box", float(F(r["box_reflected_area"][1])))
