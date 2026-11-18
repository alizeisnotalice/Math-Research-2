"""Area enclosure of the universal cover with EH depth 0.00055."""
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

def read(name):
    return json.loads(Path(__file__).with_name(name).read_bytes())

geometry = read("ah-eh-combined-550.json")
baseline = read("gibbs-1p5497-refined-18-128.json")
ah = read("ah-eh-combined-area-360.json")
eh_path = Path(__file__).with_name("eh-cap-switch-joint-550.json")
eh = json.loads(eh_path.read_bytes())["summary"]
assert geometry["passed"] and baseline["passed"] and ah["passed"]
assert sha256(eh_path.read_bytes()).hexdigest() == geometry["EH_baseline_certificate_sha256"]
assert len({F(v) for v in (geometry["sigma_degrees_exact"], baseline["sigma_degrees"],
                           ah["sigma_degrees_exact"], eh["sigma_degrees"])}) == 1
assert F(ah["AH_angle_radians_exact"]) == F(geometry["AH_cut_angle_radians_exact"])
assert F(ah["AH_pivot_exact"]) == F(geometry["AH_cut_pivot_parameter"])
assert F(eh["delta"]) == F(geometry["EH_cut_depth_radians_exact"]) == F(".00055")
lo, hi = map(F, baseline["area_interval"])
for cap in (eh["cut_gain_area_interval_exact"], ah["AH_cap_area_exact_interval"]):
    small, large = map(F, cap)
    assert 0 < small <= large
    lo -= large
    hi -= small
assert 0 < lo <= hi < F(geometry["new_area_upper_exact"])
scale = 10**39
floored, ceiled = lo*scale//1, -(-hi*scale//1)
assert F(floored,scale) < lo <= hi < F(ceiled,scale)
print(json.dumps(dict(exact=[str(lo),str(hi)], decimal=["0.%039d"%i for i in (floored,ceiled)]), indent=2))
