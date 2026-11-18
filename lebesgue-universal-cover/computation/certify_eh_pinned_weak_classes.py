"""Recheck any narrow-margin pinned AH class with the actual deepest EH chord."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

from flint import arb, ctx
from certify_ah_two_hit_budgets import geometry
from certify_ah_pinned_normals import validate


def verify(archive, ah, normal, offset, cost):
    assert ah["script_sha256"] == sha256((Path(archive)/"certify_ah_pinned_chord.py").read_bytes()).hexdigest()
    rows = ah["bins"] + ah["common_classes"]
    assert "root_conditioned_classes" not in ah
    ordinary = [arb(r["placement_margin_lower"])-cost for r in rows if arb(r["placement_margin_lower"]) > cost]
    weak = {r["reference_bin"] for r in ah["common_classes"] if not arb(r["placement_margin_lower"]) > cost}
    assert len(ordinary)+len(weak) == len(rows) and weak
    cfg = ah["parameters"]
    lines = [tuple(map(F, row)) for row in cfg["lines"]]
    far, negative, eps = map(F,cfg["physical_parameters"])
    bits = ctx.prec
    try:
        ns, _, _ = geometry(archive,lines,far=far,negative=negative,eps=eps)
        ns["rows"] = list(ns["rows"]) + [("EH-deep-chord",normal,offset)]
        out = []
        for low,high,l,r,turn,sign,orient in cfg["classes"]:
            for j in sorted(weak.intersection(range(low,high))):
                checked = validate(ns,j,F(l),F(r),turn,sign,orient,angular=12,depth=4)
                assert arb(checked["root_margin_lower"]) > 0 and arb(checked["placement_margin_lower"]) > 0
                out.append(checked)
        assert sorted(r["reference_bin"] for r in out) == sorted(weak)
    finally:
        ctx.prec = bits
    slack = min(x.lower() for x in ordinary+[arb(r["placement_margin_lower"]) for r in out])
    return slack, out
