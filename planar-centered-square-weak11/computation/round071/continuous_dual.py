#!/usr/bin/env python3
"""Exact finite rectangle arrangement checker and column oracle for R71.

No LP solver or numerical tolerance is used. Public API: check_dual(data).
Input coordinates and coefficients must be JSON integers or rational strings.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

SCHEMA = "r071-continuous-dual-v1"


def rational(x):
    if isinstance(x, bool) or not isinstance(x, (str, int)):
        raise ValueError("exact data must be integer or rational string, never float")
    return Fraction(x)


def rect(item, *, cap=False):
    if not isinstance(item.get("id"), str) or not item["id"]:
        raise ValueError("every rectangle needs a nonempty string id")
    b = tuple(map(rational, item["bounds"]))
    if len(b) != 4 or b[0] > b[1] or b[2] > b[3]:
        raise ValueError("bounds are [xmin,xmax,ymin,ymax] in ascending order")
    boundary = item.get("boundary")
    if boundary not in ("open", "closed") or cap and boundary != "closed":
        raise ValueError("explicit open/closed witness boundary; caps must be closed")
    w = rational(item["weight"])
    if w < 0:
        raise ValueError("all dual coefficients must be nonnegative")
    rhs = rational(item["rhs"]) if "rhs" in item else None
    if cap and rhs is not None and rhs < 0:
        raise ValueError("a negative cap capacity gives an inconsistent primal model")
    return {"id": item["id"], "bounds": b, "weight": w,
            "boundary": boundary, "rhs": rhs}


def atoms(edges):
    """The bounded closed hull is partitioned into points and open intervals."""
    e = sorted(set(edges))
    out = []
    for k, x in enumerate(e):
        out.append({"type": "point", "lo": x, "hi": x, "rep": x})
        if k + 1 < len(e):
            out.append({"type": "open_interval", "lo": x, "hi": e[k+1],
                        "rep": (x + e[k+1])/2})
    return out


def inside_axis(x, lo, hi, boundary):
    return lo <= x <= hi if boundary == "closed" else lo < x < hi


def bits(mask):
    while mask:
        k = (mask & -mask).bit_length() - 1
        yield k
        mask &= mask - 1


def atom_json(a):
    return {"type": a["type"], "lo": str(a["lo"]), "hi": str(a["hi"]),
            "representative": str(a["rep"])}


def check_dual(data, *, include_cells=True):
    """Return exact global result, minimum slack cell, and replayable table.

    `excluded_points` requires the mathematical input assertion nu({p})=0.
    `zero_capacity_cap_ids` may only cite supplied closed caps with exact rhs=0;
    those sets have nu-measure zero for every primal-feasible nonnegative nu.
    Without exclusions the returned certificate is pointwise on all of R^2.
    """
    start = time.perf_counter()
    if data.get("schema") != SCHEMA:
        raise ValueError("unsupported schema")
    ws = [rect(x) for x in data.get("witnesses", [])]
    cs = [rect(x, cap=True) for x in data.get("caps", [])]
    rs = ws + cs
    ids = [x["id"] for x in rs]
    if len(set(ids)) != len(ids):
        raise ValueError("rectangle ids must be globally unique")
    excluded = {tuple(map(rational, p)) for p in data.get("excluded_points", [])}
    if any(len(p) != 2 for p in excluded):
        raise ValueError("excluded_points must be pairs")
    zero_ids = set(data.get("zero_capacity_cap_ids", []))
    caps_by_id = {x["id"]: x for x in cs}
    if any(i not in caps_by_id or caps_by_id[i]["rhs"] != 0 for i in zero_ids):
        raise ValueError("zero-capacity exclusions need supplied caps with rhs exactly zero")
    zero_indices = {len(ws)+k for k, r in enumerate(cs) if r["id"] in zero_ids}
    # Closed or open degenerate witnesses are retained. Open degenerates have
    # empty support and can never produce an active witness below.
    if not ws:
        xa, ya = [], []
    else:
        hull = (min(r["bounds"][0] for r in ws), max(r["bounds"][1] for r in ws),
                min(r["bounds"][2] for r in ws), max(r["bounds"][3] for r in ws))
        xe = [v for r in rs for v in r["bounds"][:2] if hull[0] <= v <= hull[1]]
        ye = [v for r in rs for v in r["bounds"][2:] if hull[2] <= v <= hull[3]]
        xe.extend(p[0] for p in excluded if hull[0] <= p[0] <= hull[1])
        ye.extend(p[1] for p in excluded if hull[2] <= p[1] <= hull[3])
        xa, ya = atoms(xe), atoms(ye)
    masks = []
    for axis, aa in enumerate((xa, ya)):
        masks.append([sum(1 << k for k, r in enumerate(rs)
                          if inside_axis(a["rep"], r["bounds"][2*axis],
                                         r["bounds"][2*axis+1], r["boundary"]))
                      for a in aa])
    witness_mask = (1 << len(ws)) - 1
    zero_mask = sum(1 << k for k in zero_indices)
    signature_cache = {}
    cells = []
    min_cell = None
    min_slack = None
    minimum_by_dimension = {}
    negative_count = skipped = relevant = 0
    for ix, x in enumerate(xa):
        for iy, y in enumerate(ya):
            mask = masks[0][ix] & masks[1][iy]
            if not mask & witness_mask:
                continue
            relevant += 1
            if mask & zero_mask or (x["rep"], y["rep"]) in excluded:
                skipped += 1
                continue
            if mask not in signature_cache:
                wi = [k for k in bits(mask & witness_mask)]
                ci = [k for k in bits(mask & ~witness_mask)]
                slack = sum((rs[k]["weight"] for k in ci), Fraction(0)) - sum(
                    (rs[k]["weight"] for k in wi), Fraction(0))
                signature_cache[mask] = (wi, ci, slack)
            wi, ci, slack = signature_cache[mask]
            cell = {"axis_atom_indices": [ix, iy],
                    "dimension": (x["type"] != "point") + (y["type"] != "point"),
                    "cell_type": x["type"] + "_x_" + y["type"],
                    "representative": [str(x["rep"]), str(y["rep"])],
                    "active_witness_ids": [rs[k]["id"] for k in wi],
                    "active_cap_ids": [rs[k]["id"] for k in ci],
                    "slack": str(slack)}
            if include_cells:
                cells.append(cell)
            if min_slack is None or slack < min_slack:
                min_slack, min_cell = slack, cell
            dim = cell["dimension"]
            if dim not in minimum_by_dimension or slack < rational(minimum_by_dimension[dim]["slack"]):
                minimum_by_dimension[dim] = cell
            negative_count += slack < 0
    complete_rhs = all(r["rhs"] is not None for r in rs)
    constant = sum((r["weight"] * r["rhs"] for r in rs), Fraction(0)) if complete_rhs else None
    # The constant uses the caller's chosen affine witness convention. It is
    # merely arithmetic, not automatic certification of any supplied RHS.
    return {"schema": "r071-continuous-dual-result-v1", "passed": negative_count == 0,
            "scope": "pointwise" if not excluded and not zero_ids else "nu_almost_everywhere_under_declared_exclusions",
            "minimum_slack": str(min_slack) if min_slack is not None else None,
            "most_violated_cell": min_cell if negative_count else None,
            "minimum_cell": min_cell, "negative_cell_count": negative_count,
            "minimum_cells_by_dimension": {str(k):v for k,v in sorted(minimum_by_dimension.items())},
            "witness_union_cells": relevant, "checked_cells": relevant-skipped,
            "excluded_cells": skipped, "unique_checked_signatures": len(signature_cache),
            "axis_atoms": {"x": list(map(atom_json, xa)), "y": list(map(atom_json, ya))},
            "cells": cells if include_cells else None,
            "affine_constant_if_all_rhs_supplied": str(constant) if constant is not None else None,
            "exclusion_assumptions": {"points": [[str(v) for v in p] for p in sorted(excluded)],
                                     "zero_capacity_cap_ids": sorted(zero_ids)},
            "outside_witness_union": "LHS=0; RHS>=0 by nonnegative cap weights",
            "arithmetic": "Python fractions.Fraction; exact integer/rational inputs; no tolerance",
            "wall_seconds": time.perf_counter()-start}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--summary-only", action="store_true")
    args = p.parse_args()
    raw = args.input.read_bytes()
    out = check_dual(json.loads(raw), include_cells=not args.summary_only)
    out["receipt"] = {"python": sys.version, "platform": platform.platform(),
                      "input_sha256": hashlib.sha256(raw).hexdigest(),
                      "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps({k: out[k] for k in ("passed", "scope", "minimum_slack", "negative_cell_count", "checked_cells")}))


if __name__ == "__main__":
    main()
