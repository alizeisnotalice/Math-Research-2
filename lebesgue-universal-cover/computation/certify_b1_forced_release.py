"""Quantify the outward patch needed to retain a backup placement on a continuum."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb
from certify_ah_two_hit_budgets import geometry


def verify(archive):
    fold_file = archive / "ah-ideal-root-fold-certificate.json"
    obstruction = json.loads((archive / "ah-pivot-obstruction-certificate.json").read_bytes())
    fold = json.loads(fold_file.read_bytes())
    assert sha256(fold_file.read_bytes()).hexdigest() == obstruction["fold_certificate_sha256"]
    ns, _, _ = geometry(archive, [(F(".427"), F(".0048"))], far=F(".183"),
                        negative=F(".210"), eps=F(".0051"))
    A, dot, sub, perp = (ns[k] for k in ("A", "dot", "sub", "perp"))
    R = 1/ns["r3"]
    V = (A(0), R)
    normals = list(map(ns["av"], ns["d"]["n"]))
    left, right, pivot = map(F, (".42885", ".42999", ".430"))
    assert F(fold["unique_fold_exact_enclosure"][1]) < left
    assert arb(fold["discriminant_derivative"][1]) < 0
    delta = A(".018")
    steps = 114
    width = (right-left)/steps
    min_t, area = A(1), A(0)
    rows = []
    for i in range(steps):
        a, b = left+i*width, left+(i+1)*width
        vertices = ns["centered_points"](a, a)
        def coefficient(direction, indices):
            return sum((-1)**k * dot(vertices[j], direction(normals[k]))
                       for k, j in enumerate(indices))
        c = coefficient(lambda z: z, (2, 3, 0))
        d = coefficient(perp, (2, 3, 0))
        # The discriminant decreases with w and bounds the chamber maximum.
        gap = min((A(".5")-(c*c+d*d).sqrt()).lower(),
                  arb(fold["strict_negativity_margin_in_cells"]).lower())
        vertices = ns["centered_points"](a, b)
        slope = coefficient(perp, (1, 3, 0))
        identity = ns["norm"](ns["D"]) * A(pivot-b) * delta.sin()
        caps = arb(obstruction["min_other_dihedral_violation_of_actual_m3_or_m5"])
        threshold = min((gap/(6*R)).lower(), (identity/(2*R)).lower(),
                        (caps/(2*R)).lower())
        t = F(int((threshold*10**12).lower().floor().unique_fmpz())-1, 10**12)
        assert t > 0 and ((1-A(t))*slope + 3*R*A(t)
                          + A("1.5")*(A(".25")*ns["angle_unit"]).tan()) < 0
        q = ns["reflection"](vertices[2])
        dq = ns["reflection"](ns["point_jets"](ns["box"](a,b))[1][2])
        ray = sub(V,q)
        det = ns["cross"](dq,ray)
        assert ray[1]>0 and det<0
        area += (A(t)-A(t)**2/2)*(-det).lower()*A(width)
        min_t = min(min_t,A(t).lower())
        rows.append(str(t))
    index = 8
    sample = left + index*width
    mixing = F(rows[index])
    cap_slack = arb(obstruction["min_other_dihedral_violation_of_actual_m3_or_m5"])-2*R*A(mixing)
    tangent_excess = A(mixing)*arb(obstruction["min_reflected_triangle_vertex_above_true_paired_tangent"])
    assert cap_slack > A(".0009") and tangent_excess > A(".00000000577457")
    return dict(witness_parameter=str(sample), witness_mixing=str(mixing),
                witness_external_distance_lower=str(tangent_excess.lower()), pivot=str(pivot), angle_minimum=".018", parameter_range=[str(left),str(right)],
                steps=steps, permissible_weights=rows, minimum_weight=str(min_t),
                forced_external_area_lower=str(area.lower()),
                fold_sha256=obstruction["fold_certificate_sha256"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(verify(Path(__file__).resolve().parent), indent=2)+"\n")
