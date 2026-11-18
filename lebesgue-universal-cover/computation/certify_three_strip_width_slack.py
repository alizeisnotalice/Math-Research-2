"""Small width-slack strengthening of the audited exact-three-strip obstruction."""
from fractions import Fraction as F
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

import rational_trig as trig

base = Path(__file__).with_name("certify_three_strip_cap_families.py")
with contextlib.redirect_stdout(io.StringIO()):
    env = runpy.run_path(str(base))
record=env["harmonics"]
tolerance, cap_slack=env["tolerance"],env["cap_slack"]
relative_strip_half_slack=F(1,100000)

def sine_degrees(q):
    return trig.trig(trig.scale(trig.PI,F(q)/180),"sin")

def lower(x):
    return x[0]

phases={}
for k in (7,13):
    item=record[str(k)]
    relative_S_error=F(item["relative_S_minus_linear_upper"])
    actual_residue_norm_lower=(F(item["exact_base_residue_norm_lower"])
                               -F(item["residue_perturbation_absolute_upper"]))
    rho=F(1,4*(k*k-1))
    phase_allowance=3*relative_strip_half_slack/(rho*(k-1)*F(3,4)*(1-relative_S_error))
    assert phase_allowance<lower(sine_degrees(F(1,2)))
    translation_error=5*relative_strip_half_slack*tolerance
    excess=rho*actual_residue_norm_lower*lower(sine_degrees(F(3,2)))-translation_error
    assert excess>cap_slack
    phases[str(k)]={
        "relative_phase_under_relaxed_strips_upper":str(phase_allowance),
        "extra_orientation_phase_error_strictly_below_degrees":"1/2",
        "absolute_translated_cap_support_error_at_most":str(translation_error),
        "relaxed_strips_violating_support_excess_lower":str(excess),
    }
print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "exact_strip_base_source_sha256":hashlib.sha256(base.read_bytes()).hexdigest(),
    "strip_half_slack_per_unit_max_angular_deviation":str(relative_strip_half_slack),
    "each_relative_normal_deviation_radians_at_most":str(tolerance),
    "each_cap_outward_slack_at_most":str(cap_slack),
    "phase_audit":phases,
    "passed":True,
},indent=2))
