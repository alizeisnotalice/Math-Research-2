"""Bound the six-body minimum away from the five-body contact placement."""
from fractions import Fraction as F
import json
from pathlib import Path

def read(name):
    return json.loads(Path(__file__).with_name(name).read_bytes())

box = read("basin-083945-input.json")
wide = read("wide-083730-input.json")
basin = read("five-body-contact-basin.json")
fit = read("competing-fits-input.json")
assert box["center"] == wide["center"] == basin["raw_center"]
assert box["halfwidths"] == basin["raw_halfwidth"]
assert box["terms"] == wide["terms"] == fit["terms"]
assert set(map(F, wide["halfwidths"])) == {F(1, 1200)}
old = list(map(F, box["halfwidths"]))
assert old == [F(".00000025")]*2 + [F(".000055")]*9
half = F(".000377")
sqrt2_upper = F(1414214, 1000000)
assert sqrt2_upper**2 > 2
error = max((half-old[0])*sqrt2_upper, (half-old[2])*(sqrt2_upper+F(".163")))
core = F(read("basin-083945-parent.json")["lower"])
area_floor = (core-F(22,7)*error**2)/(1+4*error)
upper = F(read("competing-fits-replay.json")["five_with_s_strictly_inside"]["area_upper"])
assert area_floor-upper > F(".0000065")
print("six-body gap on the enlarged box exceeds", F(".0000065"))
