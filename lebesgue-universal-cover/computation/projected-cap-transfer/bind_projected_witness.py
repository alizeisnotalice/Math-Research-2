"""Join the independent cap proof to the actual exact-unit point clouds."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from itertools import combinations
from flint import arb, ctx


def check(comp, dest):
    assert __debug__
    source = comp / "projected-cap-transfer" / "verify.py"
    pinned = "cf28bf2f8664b3bbc427b657f844ad1edd58ca86819e7e7246f366abc3d93e88"
    assert sha256(source.read_bytes()).hexdigest() == pinned
    spec = importlib.util.spec_from_file_location("projected_cap_proof", source)
    cap = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cap)
    data = json.loads((dest / "joint-receipt.json").read_bytes())
    rho = Q(data["rho"])
    proof = cap.certify(comp, rho)
    assert proof["complete"]
    assert all(
        Q(proof[key]) == Q(data[label])
        for key, label in (
            ("point_shrink", "rho"),
            ("scale", "s"),
            ("threshold", "G"),
            ("factor_squared", "p30_exact_max_factor_squared"),
        )
    )
    factor_sq = Q(data["factor_squared"])
    lower = Q(proof["threshold"]) * factor_sq
    assert Q(data["conditional_area_lower"]) == lower <= Q(proof["area_lower"])
    raw = (dest / "joint-clouds.json").read_bytes()
    digest = sha256(raw).hexdigest()
    assert digest == data["coordinates_sha256"]
    clouds = json.loads(raw)
    assert [g["body"] for g in clouds] == ["D", "B3", "B5", "B7", "B9", "P30"]
    assert [len(g["points"]) for g in clouds] == data["counts"]
    assert all(
        Q(r["inner_scale_squared"]) >= Q(proof["scale"]) ** 2 * factor_sq
        for r in data["group_certificates"][:-1]
    )
    checks = []
    for group in clouds:
        points = group["points"]
        sq = max(
            (x - u) ** 2 + (y - v) ** 2 for (x, y), (u, v) in combinations(points, 2)
        )
        assert sq == group["denominator"] ** 2
        checks.append(len(points) * (len(points) - 1) // 2)
    assert sum(checks) == data["integer_pairs"]
    places = 24
    q = lower.numerator * 10**places // lower.denominator
    floor = f"{q//10**places}.{q%10**places:0{places}d}"
    assert Q(floor) < lower
    ctx.prec = 192
    n = max(data["counts"])
    previous_angle = arb.pi() / (n - 1)
    target_scaled = Q(data["s"]) ** 2 * factor_sq
    prev_gap = (previous_angle.cos() / (previous_angle / 2).cos()) ** 2 - arb(
        target_scaled.numerator
    ) / target_scaled.denominator
    optimal = bool(prev_gap < 0)
    if rho == Q(".000000872664") and factor_sq == Q(
        data["p30_exact_max_factor_squared"]
    ):
        assert optimal
    joined = dict(
        previous_count_disk_scale_difference=str(prev_gap),
        point_count_optimal_for_common_disk_containment=optimal,
        complete=True,
        rho=str(rho),
        area_lower=str(lower),
        strict_decimal_lower=floor,
        point_counts=data["counts"],
        integer_pairs=sum(checks),
        all_final_diameter_squared=1,
        projection_slabs=proof["projection_slabs"],
        cap_sources=proof["sources"],
        cap_script_sha256=pinned,
        point_clouds_sha256=digest,
        finite_witness_replay_sha256=sha256(
            (dest / "joint-receipt.json").read_bytes()
        ).hexdigest(),
        joiner_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (dest / "joined-proof.json").write_text(
        json.dumps(joined, sort_keys=True, indent=2) + "\n"
    )
    print(
        joined["strict_decimal_lower"],
        joined["integer_pairs"],
        joined["projection_slabs"],
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("computation", type=Path)
    ap.add_argument("cloud_directory", type=Path)
    args = ap.parse_args()
    check(args.computation, args.cloud_directory)
