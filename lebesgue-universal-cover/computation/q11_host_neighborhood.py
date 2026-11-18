"""Transfer the eleven-point support gap to every host near the free P30 pose."""
import json
from fractions import Fraction as F
from pathlib import Path
from flint import arb, ctx

root = Path(__file__).resolve().parent
host = json.loads((root/"free-p30-reoptimized-audit.json").read_bytes())
obstruction = json.loads((root/"q11-free-p30-escape-audit.json").read_bytes())
assert host["pose_sha256"] == obstruction["host_sha256"]
baseline = F(host["area_lower"])
radius = F(host["hull_radius_upper"])
ctx.prec = 160
gap = F(".0022191386")
assert all(arb(h["area_independent_padding_lower"]) > arb(gap.numerator)/gap.denominator
           for h in obstruction["hands"].values())
result = []
for displacement, target in [(F(0), F(".83751978")), (F(".000001"), F(".83751642")),
                             (F(".000002"), F(".83751307"))]:
    r = F(1,2)-displacement
    area = (baseline-F(22,7)*displacement**2)/(1+2*displacement/r)
    offset = gap-displacement
    # On the shifted host the support wall is <=radius+displacement.
    cap_denominator2 = (radius+gap)**2-r**2
    assert target > area and r > 0 and offset > 0 and cap_denominator2 > 0
    assert (target-area)**2*cap_denominator2 < r**2*offset**4
    result.append(dict(hausdorff_radius=str(displacement),augmented_area_strict_lower=str(target)))
print(json.dumps(result,indent=2))
