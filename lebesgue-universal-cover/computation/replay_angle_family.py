"""Replay the upper-bound proof's interval inequalities over one shared angle interval."""
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

base = Path(__file__).resolve().parent
if len(sys.argv) != 3:
    raise SystemExit("usage: replay_angle_family.py SIGMA_DEG_LO SIGMA_DEG_HI")
bounds = [Fraction(value) for value in sys.argv[1:]]
assert 0 < bounds[0] < bounds[1] < 30
canonical = [str(value) for value in bounds]
env = {**os.environ, "LEBESGUE_SIGMA_DEGREES": ",".join(canonical)}

checks = {
    "AH_circle_obstacle": "certify_ah_shortest",
    "AH_contact_order": "certify_ah_contact_order",
    "EH_strict_interior": "audit_eh_geometry_final",
    "five_relocations": "certify_cover_relocations",
    "global_splice": "certify_splice_connectors",
    "area": "certify_cut_area_formula",
}
results = {}
for label, module in checks.items():
    result = subprocess.run(
        [sys.executable, "-B", str(base / (module + ".py"))],
        env=env,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        raise RuntimeError(f"{label} returned {result.returncode}\n{result.stderr}")
    report = json.loads(result.stdout)
    assert report["sigma_degrees_interval"] == canonical, label
    assert report.get("passed", report.get("exact_maps_passed")) is True, label
    results[label] = report

area = results["area"]
lo, hi = (Fraction(value) for value in area["full_formula_interval"])
assert area["strict_target_passed"] == (hi < Fraction("0.8440935944"))
sources = [
    "rational_trig", "sigma_parameter", "replay_angle_family",
    "certify_sprague_fixed", "certify_fans_in_p", "certify_ah_path",
    "certify_ah_coordinates", "certify_ah_area", "certify_ah_curvature",
    "audit_eh_geometry", *checks.values(),
]
output = {
    "python": sys.version,
    "angle_interval_degrees_exact": canonical,
    "geometry_checks_passed": True,
    "area_interval_exact": [str(lo), str(hi)],
    "strict_area_target": "0.8440935944",
    "strict_area_target_passed": area["strict_target_passed"],
    "source_hashes": {
        name + ".py": hashlib.sha256((base / (name + ".py")).read_bytes()).hexdigest()
        for name in sources
    },
    "entrypoint_receipts": results,
}
print(json.dumps(output, indent=2))
