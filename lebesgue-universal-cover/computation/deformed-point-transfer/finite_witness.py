"""Build six exact unit-diameter rational sets at a requested target below G/s^2."""

from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import gzip
import importlib.util
from itertools import combinations
import json
from pathlib import Path
from flint import arb, ctx


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def sub(p, q):
    return p[0] - q[0], p[1] - q[1]


def hull(P):
    S = sorted(set(map(tuple, P)))

    def half(seq):
        H = []
        for p in seq:
            while len(H) > 1 and det(sub(H[-1], H[-2]), sub(p, H[-1])) <= 0:
                H.pop()
            H.append(p)
        return H[:-1]

    return half(S) + half(reversed(S))


def edges(H):
    i = min(range(len(H)), key=lambda k: (H[k][1], H[k][0]))
    H = H[i:] + H[:i]
    return H[0], [sub(q, p) for p, q in zip(H, H[1:] + H[:1])]


def axis_max(P):
    # Merge the CCW edge fans of H and -H. At each vertex of H-H, update squared radius.
    H = hull(P)
    assert len(H) >= 3
    p, A = edges(H)
    q, B = edges([(-x, -y) for x, y in H])
    x, y = p[0] + q[0], p[1] + q[1]
    z = (x * x + y * y, (x, y))
    i = j = 0

    def side(u):
        return int(u[1] < 0 or (u[1] == 0 and u[0] < 0))

    while i < len(A) or j < len(B):
        if i == len(A):
            cmp = 1
        elif j == len(B):
            cmp = -1
        else:
            c = det(A[i], B[j])
            cmp = side(A[i]) - side(B[j]) or ((c < 0) - (c > 0))
        if cmp <= 0:
            x += A[i][0]
            y += A[i][1]
            i += 1
        if cmp >= 0:
            x += B[j][0]
            y += B[j][1]
            j += 1
        z = max(z, (x * x + y * y, (x, y)))
    assert (x, y) == (p[0] + q[0], p[1] + q[1])
    return z


def run(comp, proofdir, dest, threshold):
    assert __debug__
    proof = json.loads((proofdir / "proof.json").read_bytes())
    P = (proofdir / "points.json").read_bytes()
    assert proof["complete"] and sha256(P).hexdigest() == proof["points_sha256"]
    s, G, area, d2 = map(
        F, [proof[k] for k in ["scale", "threshold", "area_lower", "squared_diameter"]]
    )
    assert area == G / (s * s) and d2 < s * s and 0 < threshold < area
    source = comp / "projected-cap-transfer" / "certify_sparse_witness.py"
    digest = sha256(source.read_bytes()).hexdigest()
    assert digest == "fb7624a306c8842dfcec7a0ad23995268c3527655819a271e267eacbe047c506"
    spec = importlib.util.spec_from_file_location("geom", source)
    g = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g)
    original = g.farthest
    # Exact checks against pair enumeration on a square, two polygons and interior points.
    for vertices in [
        [(0, 0), (1, 0), (0, 1), (1, 1)],
        [(0, 0), (2, 0), (3, 2), (1, 3), (1, 1)],
        json.loads(P)["deformed"]["points"],
    ]:
        assert axis_max(vertices)[0] == original(vertices)[0]
    g.farthest = axis_max
    ctx.prec = 224
    pi = arb.pi()
    g.S = max(
        10**21,
        10**10 * ((area - threshold).denominator // (area - threshold).numerator),
    )
    target = threshold / area
    ratio = g.ar(target).sqrt()
    groups = []
    rows = []
    for n in [0, 3, 5, 7, 9]:
        if n == 0:
            cs = (ratio + (ratio**2 + 8).sqrt()) / 4
            k = max(3, int((pi / (2 * cs.acos())).upper().ceil().unique_fmpz()))
            if not k % 2:
                k += 1
            step = 2
        else:
            k = max(1, int((pi / (2 * n * ratio.acos())).upper().ceil().unique_fmpz()))
            step = 1
        while True:
            group, sq, _, axis = g.curve_points(n, k, pi)
            inner = F((g.S - 2) ** 2, sq)
            if inner > target:
                break
            k += step
        assert axis_max(group["points"])[0] == group["denominator"] ** 2
        groups.append(group)
        rows.append(
            dict(
                body=group["body"],
                points=len(group["points"]),
                grid=k,
                squared_inner_scale=str(inner),
                integer_diameter_squared=sq,
                axis=axis,
            )
        )
    item = json.loads(P)["exact_unit"]
    item["body"] = "Y30"
    assert axis_max(item["points"])[0] == item["denominator"] ** 2
    groups.append(item)
    N = max(len(v["points"]) for v in groups)
    f2 = min(1 / d2, min(F(row["squared_inner_scale"]) for row in rows) / (s * s))
    assert threshold < G * f2 < area

    def disk_upper(n):
        return ((pi / n).cos() / (pi / (2 * n)).cos()) ** 2 * g.ar(area)

    # The best possible disk bound with n vertices limits this particular transfer method.
    previous = disk_upper(N - 1)
    prevodd = disk_upper(N - 2)
    dest.mkdir(parents=True, exist_ok=True)
    raw = g.packed(groups)
    gz = gzip.compress(raw, compresslevel=7, mtime=0)
    (dest / "clouds.json.gz").write_bytes(gz)
    out = dict(
        complete=True,
        threshold=str(threshold),
        factor_squared=str(f2),
        area_lower=str(G * f2),
        strict_floor=g.decimal_down(G * f2, 24),
        counts=[len(v["points"]) for v in groups],
        max_count=N,
        group_certificates=rows,
        previous_count_disk_area_upper=str(previous),
        previous_odd_disk_area_upper=str(prevodd),
        previous_count_excluded_by_disk=bool(previous < g.ar(threshold)),
        two_fewer_excluded_by_disk=bool(prevodd < g.ar(threshold)),
        grid_scale=g.S,
        coordinates_sha256=sha256(raw).hexdigest(),
        source_sha256=digest,
        deformation_sha256=proof["points_sha256"],
        verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (dest / "finite-proof.json").write_text(
        json.dumps(out, sort_keys=True, indent=2) + "\n"
    )
    print(
        out["counts"],
        out["strict_floor"],
        out["previous_count_excluded_by_disk"],
        flush=True,
    )


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    for k in ["computation", "proof", "destination"]:
        p.add_argument(k, type=Path)
    p.add_argument("threshold", type=F)
    a = p.parse_args()
    run(a.computation, a.proof, a.destination, a.threshold)
