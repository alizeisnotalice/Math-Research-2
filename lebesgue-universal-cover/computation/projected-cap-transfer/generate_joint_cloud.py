"""Generate exact-unit rational witnesses for a supplied shrink lemma and common squared factor."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import platform
import flint
from flint import arb, ctx


def load_geometry():
    path = Path(__file__).with_name("certify_sparse_witness.py")
    digest = sha256(path.read_bytes()).hexdigest()
    assert digest == "fb7624a306c8842dfcec7a0ad23995268c3527655819a271e267eacbe047c506"
    spec = importlib.util.spec_from_file_location("exact_sparse_geometry", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, digest


def main(comp, dest, rho, factor_squared):
    assert __debug__
    geo, source_hash = load_geometry()
    if not 0 <= rho < 1:
        raise ValueError("rho must lie in [0,1)")
    input_file = comp / "finite-point-leaf-input.json"
    assert (
        sha256(input_file.read_bytes()).hexdigest()
        == geo.PINS["finite-point-leaf-input.json"]
    )
    data = json.loads(input_file.read_bytes())
    p = data["points"]
    pdiam, paxis = geo.farthest(p)
    diam_sq = Q(pdiam, data["denominator"] ** 2)
    s = Q(".999999006")
    G = Q(120689264254125501, 2**57)
    max_factor_sq = 1 / (diam_sq * (1 - rho) ** 2)
    if factor_squared is None:
        factor_squared = max_factor_sq
    if not 0 < factor_squared <= max_factor_sq:
        raise ValueError(
            "requested squared factor exceeds the admissible exact P30 scale or is nonpositive"
        )
    target_sq = s * s * factor_squared
    if target_sq >= 1:
        raise ValueError(
            "the five target bodies must have a common diameter strictly below one"
        )
    slack = 1 - target_sq
    digits = max(15, len(str(10**6 * slack.denominator // slack.numerator)))
    geo.S = 10**digits
    ctx.prec = max(224, 4 * digits + 100)
    pi = arb.pi()
    goal = geo.ar(target_sq).sqrt()
    groups = []
    rows = []
    for n in (0, 3, 5, 7, 9):
        if n == 0:
            cosine = (goal + (goal * goal + 8).sqrt()) / 4
            bound = pi / (2 * cosine.acos())
            k = max(3, int(bound.upper().ceil().unique_fmpz()))
            if k % 2 == 0:
                k += 1
            step = 2
        else:
            bound = pi / (2 * n * goal.acos())
            k = max(1, int(bound.upper().ceil().unique_fmpz()))
            step = 1
        # The count above is a proposal. Only this exact containment test authorizes output.
        while True:
            cloud, sq, _, axis = geo.curve_points(n, k, pi)
            scale_sq = Q((geo.S - 2) ** 2, sq)
            if scale_sq >= target_sq:
                break
            k += step
        groups.append(cloud)
        rows.append(
            dict(
                body=cloud["body"],
                count=len(cloud["points"]),
                grid=k,
                inner_scale_squared=str(scale_sq),
                diametral_integer_axis=axis,
                integer_diameter_squared=sq,
                final_squared_diameter="1",
            )
        )
    groups.append(geo.exact_unit_cloud("P30", p, pdiam, paxis))
    rows.append(
        dict(
            body="P30",
            count=len(p),
            original_squared_diameter=str(diam_sq),
            normalized_squared_scale=str(1 / diam_sq),
            diametral_integer_axis=paxis,
            final_squared_diameter="1",
        )
    )
    lower = G * factor_squared
    all_scales = min(Q(v["inner_scale_squared"]) for v in rows[:-1])
    assert all_scales >= target_sq and factor_squared <= max_factor_sq
    max_rho = 1 - geo.ar(s) / (geo.ar(diam_sq * all_scales)).sqrt()
    raw = geo.packed(groups)
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "joint-clouds.json").write_bytes(raw)
    out = dict(
        schema="R221-parametric-joint-cloud-v1",
        hypothesis="Every independent placement of the five s-scaled curves and the (1-rho)-scaled P30 has hull area at least G; the generator checks only the finite-witness consequence.",
        rho=str(rho),
        s=str(s),
        G=str(G),
        factor_squared=str(factor_squared),
        p30_exact_max_factor_squared=str(max_factor_sq),
        target_curved_scale_squared=str(target_sq),
        conditional_area_lower=str(lower),
        conditional_area_floor=geo.decimal_down(lower, 21),
        same_clouds_P30_limited_max_rho=str(max_rho),
        counts=[len(g["points"]) for g in groups],
        total=sum(len(g["points"]) for g in groups),
        integer_pairs=sum(
            len(g["points"]) * (len(g["points"]) - 1) // 2 for g in groups
        ),
        group_certificates=rows,
        grid_scale=geo.S,
        precision=ctx.prec,
        arbitrary_rigid_motion=True,
        coordinates_sha256=sha256(raw).hexdigest(),
        generator_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        geometry_sha256=source_hash,
        p30_input_sha256=geo.PINS["finite-point-leaf-input.json"],
        python=platform.python_version(),
        python_flint=flint.__version__,
    )
    (dest / "joint-receipt.json").write_text(
        json.dumps(out, sort_keys=True, indent=2) + "\n"
    )
    print(
        out["counts"],
        out["conditional_area_floor"],
        out["same_clouds_P30_limited_max_rho"],
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("computation", type=Path)
    ap.add_argument("destination", type=Path)
    ap.add_argument("--rho", type=Q, required=True)
    ap.add_argument("--factor-squared", type=Q)
    args = ap.parse_args()
    main(args.computation, args.destination, args.rho, args.factor_squared)
