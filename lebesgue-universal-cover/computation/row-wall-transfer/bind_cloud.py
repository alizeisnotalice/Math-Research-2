"""Join the area-transfer proof to the generated coordinates and audit every final pair."""

from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path


def check(comp, proof_path, dest):
    assert __debug__
    proof = json.loads(proof_path.read_bytes())
    receipt = json.loads((dest / "joint-receipt.json").read_bytes())
    raw = (dest / "joint-clouds.json").read_bytes()
    assert proof["complete"] and proof["joint_leaf_row_pairs"] == 3273 * 2493
    assert sha256(raw).hexdigest() == receipt["coordinates_sha256"]
    for k in ["generate_joint_cloud.py", "certify_sparse_witness.py"]:
        key = "generator_sha256" if k.startswith("generate") else "geometry_sha256"
        assert (
            sha256((comp / "projected-cap-transfer" / k).read_bytes()).hexdigest()
            == receipt[key]
        )
    assert (
        proof["source_hashes"]["finite-point-leaf-input.json"]
        == receipt["p30_input_sha256"]
    )
    for a, b in [
        ("scale", "s"),
        ("point_shrink", "rho"),
        ("threshold", "G"),
        ("factor_squared", "factor_squared"),
        ("lower", "conditional_area_lower"),
    ]:
        assert F(proof[a]) == F(receipt[b])
    groups = json.loads(raw)
    assert [g["body"] for g in groups] == ["D", "B3", "B5", "B7", "B9", "P30"]
    assert [len(g["points"]) for g in groups] == receipt["counts"]
    assert F(receipt["p30_exact_max_factor_squared"]) == F(proof["factor_squared"])
    factor = F(proof["factor_squared"])
    assert all(
        F(g["inner_scale_squared"]) >= F(proof["scale"]) ** 2 * factor
        for g in receipt["group_certificates"][:-1]
    )
    assert F(proof["point_diameter_squared"]) == F(
        receipt["group_certificates"][-1]["original_squared_diameter"]
    )
    pair_count = 0
    for group in groups:
        points = group["points"]
        sq = max(
            (x - u) ** 2 + (y - v) ** 2 for (x, y), (u, v) in combinations(points, 2)
        )
        assert sq == group["denominator"] ** 2
        pair_count += len(points) * (len(points) - 1) // 2
    assert pair_count == receipt["integer_pairs"]
    lower = F(proof["lower"])
    assert lower == factor * F(proof["threshold"]) > F(proof["strict_floor"])
    result = dict(
        complete=True,
        area_lower=str(lower),
        floor=proof["strict_floor"],
        diameters_squared=[1] * 6,
        point_counts=receipt["counts"],
        integer_pairs=pair_count,
        max_count=max(receipt["counts"]),
        cloud_sha256=sha256(raw).hexdigest(),
        cloud_receipt_sha256=sha256(
            (dest / "joint-receipt.json").read_bytes()
        ).hexdigest(),
        row_proof_sha256=sha256(proof_path.read_bytes()).hexdigest(),
        joiner_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (dest / "joined-proof.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n"
    )
    print(result["floor"], pair_count, result["max_count"])


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("computation", type=Path)
    p.add_argument("proof", type=Path)
    p.add_argument("cloud_directory", type=Path)
    a = p.parse_args()
    check(a.computation, a.proof, a.cloud_directory)
