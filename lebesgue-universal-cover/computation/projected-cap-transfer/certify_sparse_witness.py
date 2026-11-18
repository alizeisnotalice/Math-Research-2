"""Same-pose R221 finite witnesses with a strict, exact integer diameter test for each cloud."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import platform
import flint
from flint import arb, ctx

S = 10**15
PINS = {
    "finite-point-leaf-input.json": "b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc",
    "six-input-083745/points/normal-gap-083745-audit.json": "2b76fccefe1e5d56c0752ed077957af2fb0067e56c494a67fbc8ba836fab7115",
    "six-input-083745/six-input-transfer.md": "f31a6e1dc56ce41ceac1e51800eb01b6657b91ad21abc5e9728169fdebe8a2b7",
}


def ar(x):
    r = Q(x)
    return arb(r.numerator) / r.denominator


def packed(obj):
    return (json.dumps(obj, sort_keys=True, separators=(",", ":")) + "\n").encode()


def nearest(x):
    z = (S * x + ar(Q(1, 2))).floor().unique_fmpz()
    assert z is not None
    assert abs(S * x - int(z)) < ar(Q(1, 2))
    return int(z)


def unit(angle):
    sn, cs = angle.sin_cos()
    return cs, sn


def farthest(points):
    maximum = -1
    axis = None
    for (x, y), (u, v) in combinations(points, 2):
        dx, dy = x - u, y - v
        sq = dx * dx + dy * dy
        if sq > maximum:
            maximum = sq
            axis = (dx, dy)
    assert maximum > 0
    return maximum, axis


def exact_unit_cloud(body, points, diameter, axis):
    # Division by the complex integer a+ib normalizes diameter without irrational coordinates.
    a, b = axis
    assert a * a + b * b == diameter
    coords = [[a * x + b * y, a * y - b * x] for x, y in points]
    return dict(body=body, points=coords, denominator=diameter)


def curve_points(n, k, pi):
    points = []
    if n == 0:
        assert k % 2 == 1
        R = 1 / (2 * (pi / k).cos())
        for j in range(k):
            a, b = unit(2 * pi * j / k)
            points.append([nearest(R * a), nearest(R * b)])
    else:
        R = 1 / (2 * (pi / (2 * n)).cos())
        assert 1 - R > ar(Q(2, 5))
        T = 1 / (pi / (2 * n * k)).cos()
        for i in range(n):
            ai = 2 * pi * i / n
            a, b = unit(ai)
            x, y = R * a, R * b
            points.append([nearest(x), nearest(y)])
            for j in range(k):
                a, b = unit(
                    ai + pi * (1 - ar(Q(1, 2 * n)) + ar(Q(2 * j + 1, 2 * n * k)))
                )
                points.append([nearest(x + T * a), nearest(y + T * b)])
    assert len(points) == len(set(map(tuple, points))) == (k if n == 0 else n * (k + 1))
    sq, axis = farthest(points)
    den = isqrt(sq) + 1
    assert sq < den**2
    return exact_unit_cloud("D" if n == 0 else f"B{n}", points, sq, axis), sq, den, axis


def decimal_down(z, places=18):
    scale = 10**places
    digit = z.numerator * scale // z.denominator
    return f"{digit//scale}.{digit%scale:0{places}d}"


def main(comp, dest):
    assert __debug__
    ctx.prec = 224
    pi = arb.pi()
    for name, dig in PINS.items():
        assert sha256((comp / name).read_bytes()).hexdigest() == dig
    proof = json.loads(
        (comp / "six-input-083745/points/normal-gap-083745-audit.json").read_bytes()
    )
    G, s, rho = map(
        Q, (proof["underlying_goal"], proof["scale"], proof["point_shrink"])
    )
    assert (G, s, rho) == (
        Q(120689264254125501, 2**57),
        Q(".999999006"),
        Q(".000000867442"),
    )
    P = json.loads((comp / "finite-point-leaf-input.json").read_bytes())
    p0 = P["points"]
    den0 = P["denominator"]
    pdiam, paxis = farthest(p0)
    psq = Q(pdiam, den0**2)
    t = Q(isqrt(10**36 * psq.denominator // psq.numerator), 10**18)
    pscale = t / (1 - rho)
    assert pscale == Q(proof["factor"]) and t * t * psq < 1
    sixth = exact_unit_cloud("P30", p0, pdiam, paxis)
    dest.mkdir(parents=True, exist_ok=True)
    out = {}
    for size, goal in (
        (1931, ".8374500023"),
        (2919, ".83745093731"),
        (5373, ".83745145012"),
        (5425, ".83745145414617248855"),
    ):
        groups = []
        records = []
        squared_factors = [1 / (psq * (1 - rho) ** 2)]
        # To keep the polygon groups small, choose the least arc grid (in the analytic model)
        # whose containment factor attains the limiting disk factor or the P30 factor.
        disk_model = (pi / size).cos() / (pi / (2 * size)).cos()
        limit = min(disk_model, ar(s * pscale))
        for n in (0, 3, 5, 7, 9):
            if n:
                k_bound = (pi / (2 * n * limit.acos())).ceil().unique_fmpz()
                assert k_bound is not None
                k = int(k_bound)
                assert (pi / (2 * n * k)).cos() > limit and (
                    pi / (2 * n * (k - 1))
                ).cos() < limit
            else:
                k = size
            cloud, sq, den, axis = curve_points(n, k, pi)
            scale_squared = Q((S - 2) ** 2, sq)
            squared_factors.append(scale_squared / s**2)
            groups.append(cloud)
            records.append(
                dict(
                    body=cloud["body"],
                    count=len(cloud["points"]),
                    grid=k,
                    squared_diameter_numerator=sq,
                    inner_scale_denominator=den,
                    unit_cloud_denominator=cloud["denominator"],
                    diametral_integer_axis=axis,
                    guaranteed_body_scale_squared=str(scale_squared),
                    exact_squared_diameter="1",
                )
            )
        groups.append(sixth)
        f_squared = min(squared_factors)
        lower = G * f_squared
        assert lower > Q(goal)
        assert max(len(g["points"]) for g in groups) == size
        assert 1 / psq >= f_squared * (1 - rho) ** 2
        assert f_squared * s**2 <= min(
            Q(r["guaranteed_body_scale_squared"]) for r in records
        )
        coord = packed(groups)
        path = dest / f"clouds-{size}.json"
        path.write_bytes(coord)
        out[str(size)] = dict(
            counts=[len(g["points"]) for g in groups],
            point_total=sum(len(g["points"]) for g in groups),
            G=str(G),
            joint_scale_factor_squared=str(f_squared),
            area_lower=str(lower),
            area_floor=decimal_down(lower),
            strict_target=goal,
            curve_records=records,
            coordinates_sha256=sha256(coord).hexdigest(),
            p30_scale=str(t),
            p30_inner_comparison_sq=str(t * t * psq),
            p30_diametral_integer_axis=paxis,
            p30_final_squared_diameter="1",
            unordered_pairs=sum(
                len(g["points"]) * (len(g["points"]) - 1) // 2 for g in groups
            ),
        )
        print(size, out[str(size)]["counts"], decimal_down(lower), flush=True)
    cap = {}
    for size, target in (
        (1929, Q(".83745")),
        (5371, Q(".83745145")),
        (5424, Q(proof["area_lower"])),
    ):
        bound = ar(G) * ((pi / size).cos() / (pi / (2 * size)).cos() / ar(s)) ** 2
        assert bound < ar(target)
        cap[str(size)] = dict(
            disk_containment_area_ceiling=str(bound.upper()), target=str(target)
        )
    # An exact-unit P30 has a slightly larger scale than the 18-digit factor in the old certificate.
    assert Q(out["5425"]["area_lower"]) > Q(proof["area_lower"])
    receipt = dict(
        tiers=out,
        disk_containment_method_ceilings=cap,
        coordinate_scale=S,
        upstream_sha256=PINS,
        script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        python=platform.python_version(),
        python_flint=flint.__version__,
        precision=ctx.prec,
    )
    (dest / "sparse-receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n"
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("computation", type=Path)
    ap.add_argument("destination", type=Path)
    args = ap.parse_args()
    main(args.computation, args.destination)
