#!/usr/bin/env python3
"""Read-only integrity audit for the returned n=24 tau=2 terminal trace.

The raw trace is streamed: its group identities and distance histograms are
checked as integers.  The necessary-cone arithmetic uses Fraction exactly on
the stored decimal aggregate; it is not a sampling-accuracy statement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path


BASE = 0xA022000000000000


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def finite(x) -> bool:
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def main(raw: Path, sidecar: Path, source: Path, aggregate: Path, output: Path) -> None:
    side = json.loads(sidecar.read_text())
    assert sha(raw) == side["raw_sha256"]
    assert sha(source) == side["source_sha256"]
    meta = end = None
    records = {}
    total_hist = total_w = total_l2 = 0
    global_hist = None
    with raw.open() as f:
        for lineno, line in enumerate(f, 1):
            x = json.loads(line)
            if lineno == 1:
                assert x.get("kind") == "metadata"
                meta = x
                n, P, R = x["n"], x["particles"], x["replicates"]
                N, groups, bins = n * n, (P + 255) // 256, 2 * n - 1
                global_hist = [0] * bins
                continue
            if x.get("kind") == "run_end":
                assert end is None and finite(x.get("wall_seconds"))
                end = x
                continue
            assert x.get("kind") == "replicate" and x.get("tau") == 2
            r = x["replicate"]
            assert r not in records and 0 <= r < R
            domain = int(meta["domain"], 16)
            want = (BASE + domain + n * 1_000_000 + 2 * 10_000 + r) & ((1 << 64) - 1)
            assert x["seed_hex"].lower().zfill(16) == f"{want:016x}"
            assert not x["dead"]
            for key in ("Z", "logZ", "U_over_Z", "D", "F", "G", "F_minus_2G_minus_4D"):
                assert finite(x.get(key))
            keys = ("tau2_lift_L_groups", "tau2_lift_L_distance2_groups",
                    "tau2_lift_L_F0_groups", "tau2_lift_L_G0_groups",
                    "tau2_lift_group_LF_minus_2LG_minus_4LD",
                    "tau2_distance_L_hist_groups", "tau2_distance_L2_groups")
            assert all(isinstance(x.get(k), list) and len(x[k]) == groups for k in keys)
            rw = rh2 = 0
            for g in range(groups):
                cap = min(256, P - 256 * g)
                w = x[keys[0]][g]; dd = x[keys[1]][g]; ff = x[keys[2]][g]; gg = x[keys[3]][g]
                assert all(type(q) is int for q in (w, dd, ff, gg, x[keys[4]][g]))
                assert cap <= w <= cap * (N - 3)
                assert ff - 2 * gg == 4 * dd == 4 * int(dd)
                assert x[keys[4]][g] == 0
                row = x[keys[5]][g]
                q2 = x[keys[6]][g]
                assert isinstance(row, list) and len(row) == bins and all(type(q) is int for q in row)
                assert row[0] == row[1] == 0 and sum(row) == w
                assert type(q2) is int and w <= q2 <= (N - 3) * w
                for d, q in enumerate(row):
                    assert 0 <= q <= cap * (N - 3)
                    global_hist[d] += q
                rw += w; rh2 += q2
            total_w += rw
            total_l2 += rh2
            records[r] = {"seed_hex": x["seed_hex"], "sum_L": rw, "sum_L2": rh2}
    assert meta is not None and end is not None
    assert sorted(records) == list(range(meta["replicates"]))
    assert side["expected_metadata"] == {k: meta[k] for k in ("n", "particles", "replicates", "tau_lo", "tau_hi", "domain", "method", "move_steps")}
    a = json.loads(aggregate.read_text())
    n = meta["n"]
    h = list(map(Fraction, a["distance_L_weighted_probabilities"]))
    assert len(h) == 2 * n - 1 and h[0] == h[1] == 0 and all(q >= 0 for q in h)
    h = [q / sum(h) for q in h]
    b_over_a = Fraction(a["tau2_normalizer_mean"]) / Fraction(a["lift_full_normalizer_mean"])
    rhs = sum((Fraction(d, 2 * (d - 1)) * h[d] for d in range(2, len(h))), Fraction())
    slack = rhs - b_over_a
    cones = {}
    for r in (2, 6, 12, 18):
        tail = sum(h[r + 1:], Fraction())
        upper = Fraction(2 * r, r - 1) * (1 - b_over_a)
        cones[str(r)] = {"tail": str(tail), "necessary_upper": str(upper),
                          "slack": str(upper - tail), "slack_positive": upper - tail >= 0}
    out = {
        "scope": "Raw-integer and replay-integrity audit; cone values are exact Fraction calculations on stored decimal plug-ins, not accuracy, ESS, or mixing certificates.",
        "integrity_status": "passed",
        "cone_plugin_status": "satisfies_linear_cone" if slack >= 0 else "violates_linear_cone",
        "raw_sha256": sha(raw), "source_sha256": sha(source), "metadata": meta,
        "run_end": end, "replicate_count": len(records), "replicates": records,
        "integer_checks": {"all_group_F_minus_2G_equals_4D": True,
                           "all_group_residuals_zero": True,
                           "all_distance_rows_sum_to_L": True,
                           "total_L": total_w, "total_L2": total_l2,
                           "global_distance_L_histogram": global_hist,
                           "global_distance_histogram_sum": sum(global_hist)},
        "cone_fraction": {"B_over_A_hat": str(b_over_a), "linear_rhs": str(rhs),
                          "linear_slack": str(slack), "linear_slack_positive": slack >= 0,
                          "thresholds": cones},
        "ess_warning": {"tau2_weight": a["tau2_weight"], "liftZ_weight": a["liftZ_weight"],
                        "interpretation": "Near-one ESS/maxshare near one report severe replicate-weight concentration. They invalidate ordinary replicate-SE/mixing interpretations; no conditional-law accuracy claim follows."},
    }
    output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"integrity_status": "passed", "cone_plugin_status": out["cone_plugin_status"], "raw_sha256": out["raw_sha256"],
                      "cone_slack_positive": out["cone_fraction"]["linear_slack_positive"]}))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("raw", type=Path); p.add_argument("sidecar", type=Path)
    p.add_argument("source", type=Path); p.add_argument("aggregate", type=Path)
    p.add_argument("output", type=Path)
    q = p.parse_args(); main(q.raw, q.sidecar, q.source, q.aggregate, q.output)
