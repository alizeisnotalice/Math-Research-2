#!/usr/bin/env python3
"""Final addendum to the fixed-angle Definition-9 EH interval audit.

This deliberately leaves ``audit_eh_geometry.py`` frozen.  It reruns that
certificate and proves the extra slab, boundary-extremum, and strict-interior
statements used by the final write-up.
"""
from __future__ import annotations

import contextlib
from fractions import Fraction as F
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

base = Path(__file__).parent
with contextlib.redirect_stdout(io.StringIO()):
    d = runpy.run_path(str(base / "audit_eh_geometry.py"))

I = d["I"]
add, sub, mul, dot, norm2, perp = (d[k] for k in ("add", "sub", "mul", "dot", "norm2", "perp"))
cross = d["cross"]
T, U, J, S = (d[k] for k in ("T", "U", "J", "S"))
B3, C3, pal = (d[k] for k in ("B3", "C3", "pal"))


def out(x):
    return x.out()


# Boundary order for the removed Jordan region is T--U (chord), U--J (C3
# CCW arc), J--T (B3 CCW arc).  On either CCW circle dx/dtheta=-(x-c)_y.
# Both radial y values are strictly negative, and the short-arc dot tests in
# the frozen certificate make the derivative sign constant along its arc.
rC_U, rC_J = sub(U, C3), sub(J, C3)
rB_J, rB_T = sub(J, B3), sub(T, B3)
assert rC_U[1].hi < 0 and rC_J[1].hi < 0
assert rB_J[1].hi < 0 and rB_T[1].hi < 0
assert cross(rC_U, rC_J).lo > 0 and dot(rC_U, rC_J).lo > 0
assert cross(rB_J, rB_T).lo > 0 and dot(rB_J, rB_T).lo > 0

# Thus x is strictly increasing on U->J and J->T.  It is affine on TU.
# These interval inequalities certify the whole closed EH region's x slab:
# extrema of an affine functional on a compact planar Jordan region are on
# its boundary, and the boundary extrema are U and T by the preceding signs.
assert U[0].lo > -I(F(1, 2)).hi
assert T[0].hi < -I(F(499, 1000)).lo
assert J[0].hi < -I(F(499, 1000)).lo

# Exact Pál-line equality at U is analytic, not a decimal inference:
# U=proj_pal(C3) - sqrt(...) perp(pal), so pal.U=1/2 because |pal|=1 and
# pal.perp(pal)=0.  The interval calculation is only a consistency enclosure.
pal_U = dot(pal, U)
assert pal_U.lo < I(F(1, 2)).lo < pal_U.hi

# The frozen certificate proves pal.x>=1/2 on the complete boundary, with
# equality only at U (strictly positive arc derivatives and chord endpoint
# values).  The affine-extremum principle therefore promotes it to the closed
# Jordan region and its strict version to the open EH interior.

# S is outside the closure already from the x slab: S_x=1/2 exactly, whereas
# all of closure(EH) has x<-0.499.  For f(x)=|x-S|^2, the only possible
# interior local minimizer is S (grad f=2(x-S)).  Boundary f>=1, with equality
# only at T, was interval-certified in the frozen script.  If an interior
# point had f<=1, compactness would give an interior global minimizer; it
# would have to be S, contradiction.  Hence f>1 in EH interior.
assert S[0].lo < I(F(1, 2)).lo < S[0].hi
assert S[0].lo - T[0].hi > I(F(999, 1000)).lo

result = {
    "python": sys.version,
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "frozen_base_script_sha256": hashlib.sha256((base / "audit_eh_geometry.py").read_bytes()).hexdigest(),
    "sigma_degrees_exact": "1.5494",
    "endpoint_x_intervals": {"U": [out(U[0])], "J": [out(J[0])], "T": [out(T[0])]},
    "arc_x_monotonicity": {
        "U_to_J_C3_radial_y": [out(rC_U[1]), out(rC_J[1])],
        "J_to_T_B3_radial_y": [out(rB_J[1]), out(rB_T[1])],
        "conclusion": "x strictly increases on both CCW boundary arcs; it is affine on TU.",
    },
    "whole_closed_EH_x_slab": {
        "lower": "-1/2 < x",
        "upper": "x < -499/1000",
        "justification": "boundary monotonicity plus the affine-extremum principle on the compact Jordan closure",
    },
    "pal_line_equality": {
        "normal": "u(7pi/6)",
        "analytic_identity": "pal.U=1/2 exactly from U=proj_pal(C3)-sqrt(1-d^2)perp(pal)",
        "interval_consistency": out(pal_U),
        "closed_region": "pal.x>=1/2",
        "interior": "pal.x>1/2",
    },
    "S_disk_strict_interior": {
        "S_x_exact": "1/2",
        "S_outside_closed_EH": True,
        "boundary_fact_from_frozen_certificate": "|x-S|^2>=1 on boundary, with equality only at T",
        "critical_point_argument": "the only interior local minimizer of |x-S|^2 is S; since S is outside closure(EH), compactness implies |x-S|^2>1 throughout EH interior",
    },
    "passed": True,
    "scope": "fixed-angle geometric addendum only; does not prove the global universal-cover theorem",
}
print(json.dumps(result, indent=2))
