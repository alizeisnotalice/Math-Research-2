"""Combine three certified area enclosures using rational endpoints."""
from fractions import Fraction
from pathlib import Path
import json

def read(name):
    return json.loads(Path(__file__).with_name(name).read_bytes())

base = read("gibbs-1p5497-refined-18-128.json")
eh = read("eh-cap-switch-360.json")["summary"]
ah = read("ah-eh-combined-area-360.json")
lower, upper = map(Fraction, base["area_interval"])
for cap in (eh["cut_gain_area_interval_exact"], ah["AH_cap_area_exact_interval"]):
    small, large = map(Fraction, cap)
    assert 0 < small <= large
    lower -= large
    upper -= small
old_lo, old_hi = map(Fraction, ah["combined_candidate_area_exact_interval"])
assert old_lo <= lower <= upper <= old_hi
scale = 10**39
print("0.%039d < area < 0.%039d" % (lower*scale//1, -((-upper*scale)//1)))
