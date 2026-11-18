"""Verify the surviving G--K arc used in the quantitative mixed-area slack lemma."""
import json
from fractions import Fraction
from pathlib import Path

import audit_optimized_eh_chord as g
from certify_eh_third_facet import geometry

def main():
    comp = Path(__file__).resolve().parent
    g.SIGMA = g.arb("1.5497") * g.PI / 180
    n, m, *_ = g.setup()
    center = g.line_pair(n[3], m[2])
    first = g.circle_line(center, n[5], 1)
    last = g.minus(center, g.scale(g.dot(center, m[5]) - g.HALF, m[5]))
    a, b = g.minus(first, center), g.minus(last, center)
    # K is the projection of F3 from the opposite unit support line.
    assert all(v[0] > 0 and v[1] < 0 for v in (a, b))
    assert g.cross(a, b) > 0 and g.dot(a, b) > 0
    angle = (g.cross(a, b) / g.dot(a, b)).atan()
    assert angle > g.arb("0.014")
    assert first[0] > g.arb("0.37") and last[1] < g.arb("-0.34")
    old = json.loads((comp / "eh-two-line-cover-certificate.json").read_text())
    _, _, corners, _, _ = geometry(old)
    assert all(v[0] < 0 for v in corners)
    lower = Fraction(14, 1000)**3 / (56 * (512 + 1)**2)
    return {"sigma": "1.5497", "arc": "G--K about F3", "radius": 1,
            "angle_interval": str(angle), "angle_lower": "7/500",
            "arc_x_lower": "37/100", "arc_y_upper": "-34/100",
            "maximum_vertices_example": 512, "mixed_area_gap_strict_lower": str(lower)}

if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
