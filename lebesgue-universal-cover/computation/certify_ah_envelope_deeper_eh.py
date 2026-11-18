"""Carry the certified AH two-line placement to EH depth 0.0005521."""
from fractions import Fraction as F
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

from flint import arb, ctx

ctx.prec = 160
here = Path(__file__).resolve().parent
old_depth, new_depth = "0.00055", "0.0005521"
sigma = "1.5497"
with contextlib.redirect_stdout(io.StringIO()):
    d = runpy.run_path(str(here / "ah_root_geometry_at_angle.py"),
                       init_globals={"sigma_degrees": sigma, "eh_depth": old_depth})["d"]


def ball(q):
    if isinstance(q, d["I"]):
        return ball(q.lo).union(ball(q.hi))
    q = F(q)
    return arb(q.numerator)/q.denominator


def read(name):
    raw = (here / name).read_bytes()
    hashes[name] = hashlib.sha256(raw).hexdigest()
    return json.loads(raw)


add, sub, mul, perp, dot = (d[k] for k in ("add","sub","mul","perp","dot"))
center, U, T = (tuple(map(ball, d[k])) for k in ("c3","u","t_eh"))
assert dot(T,T) < ball(F(1,3))


def outward_normal(depth):
    x = ball(depth)
    radial = sub(U,center)
    moved = add(center,add(mul(x.cos(),radial),mul(-x.sin(),perp(radial))))
    chord = sub(T,moved)
    return mul(-1/dot(chord,chord).sqrt(),perp(chord))


hashes = {}
old_cap = read("ah-two-line-cap.json")
assert F(old_cap["sigma_degrees_exact"]) == F(sigma)
assert F(old_cap["EH_depth_radians_exact"]) == F(old_depth)
difference = sub(outward_normal(new_depth),outward_normal(old_depth))
support_change = (dot(difference,difference)*ball(F(4,3))).sqrt()
margin = arb(old_cap["uniform_actual_placement_slack_lower"])
assert margin - support_change > ball("0.000125")

old_gibbs = read("gibbs-1p5497-refined-18-128.json")
deep_eh = read("eh-cap-switch-joint-552p1-revisit-ah.json")["summary"]
ah = read("ah-two-line-area.json")
assert F(deep_eh["sigma_degrees"]) == F(old_gibbs["sigma_degrees"]) == F(sigma)
assert F(ah["sigma_degrees_exact"]) == F(sigma)
assert ah["AH_cut_parameters_exact"] == old_cap["AH_cut_parameters_exact"]
assert F(deep_eh["delta"]) == F(new_depth) and deep_eh["AH_relocation_after_root"]
begin, end = map(F,old_gibbs["area_interval"])
ah_lo, ah_hi = map(F,ah["AH_cap_area_exact_interval"])
eh_lo, eh_hi = map(F,deep_eh["cut_gain_area_interval_exact"])
lo, hi = begin-ah_hi-eh_hi, end-ah_lo-eh_lo
lower = "0.844093401110543697839148056785199639105"
upper = "0.844093401110543697839148056785199639190"
assert F(lower) < lo <= hi < F(upper)

print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "input_sha256":hashes,
    "old_EH_depth_exact":old_depth, "new_EH_depth_exact":new_depth,
    "uniform_EH_support_change_upper":str(support_change.upper()),
    "new_uniform_actual_placement_slack_lower":str((margin-support_change).lower()),
    "combined_candidate_area_exact_interval":[str(lo),str(hi)],
    "combined_candidate_area_decimal_interval":[lower,upper],
    "passed":True,
},indent=2))
