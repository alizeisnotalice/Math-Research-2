"""Convert the freshly recomputed all-fan margins into a rational area bound."""
from fractions import Fraction as Q
import json
import math
from pathlib import Path
import platform

source = Path(__file__).with_name("046-full-replay.json")
report = json.loads(source.read_text())
tree = report["tree"]
parts = report["parts"]
assert report["complete"] and tree["ok"]
assert [row["part"] for row in parts] == list(range(tree["header"]["nparts"]))
assert all(row["ok"] and row["fails"] == 0 for row in parts)
assert tree["fan_leaves"] == tree["leaves"] == sum(row["leaves"] for row in parts)
assert tree["multi_leaves"] == tree["far_leaves"] == 0
base = Q(str(tree["header"]["tau"]))
margin = min(Q(math.nextafter(row["min_margin"], -math.inf)) for row in parts)
target = Q("0.83600000162452")
assert base == Q("0.836") and base + margin > target

# For an escaped translation component, the separation exceeds the half width.
# Use inscribed disks, the two exterior half disks and the connecting trapezoid.
header = tree["header"]
bounds = []
for n, x in zip(header["ns"], [0, 3, 6]):
    w = min(-Q(header["lo"][x]), Q(header["hi"][x]),
            -Q(header["lo"][x+1]), Q(header["hi"][x+1]))
    angle_upper = Q(22, 7) / (2 * n)
    cosine_lower = 1 - angle_upper**2 / 2
    r = 1 - 1 / (2 * cosine_lower)
    lower = Q(314, 100) * (Q(1, 4) + r*r) / 2 + w*(Q(1, 2) + r)
    assert w > 0 and r > 0 and lower > target
    bounds.append(str(lower))
print(json.dumps({"python": platform.python_version(), "architecture": platform.machine(),
                  "leaf_margin_lower": str(margin), "exterior_lower": bounds,
                  "target": str(target), "passed": True}, indent=2))
