#!/usr/bin/env python3
"""Exact deterministic discrepancy certificate from a failed finite cone plug-in."""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path


def rational(x: Fraction) -> dict[str, str | float]:
    return {"fraction": str(x), "decimal": float(x)}


def main(summary: Path, output: Path) -> None:
    x = json.loads(summary.read_text())
    n = int(x["n"])
    h_raw = list(map(Fraction, x["distance_L_weighted_probabilities"]))
    assert len(h_raw) == 2 * n - 1 and h_raw[0] == h_raw[1] == 0 and all(z >= 0 for z in h_raw)
    q = [z / sum(h_raw) for z in h_raw]
    assert sum(q) == 1
    b = Fraction(x["tau2_normalizer_mean"]) / Fraction(x["lift_full_normalizer_mean"])
    g = [Fraction(0) if d < 2 else Fraction(d, 2 * (d - 1)) for d in range(2 * n - 1)]
    eg = sum((q[d] * g[d] for d in range(2, len(q))), Fraction())
    delta = b - eg
    osc = max(g[2:]) - min(g[2:])
    assert osc == Fraction(22, 45), osc
    assert delta > 0
    lower = delta / (1 + osc)
    threshold = Fraction(11, 50)
    out = {
        "scope": "Deterministic plug-in discrepancy certificate. Conditional on the exact necessary cone B/A <= E_q[g], it proves at least one of the stored B/A plug-in error and TV(q_hat,q_true) is at least lower_bound. It assumes neither a CI nor mixing and makes no bug attribution.",
        "summary_sha256": __import__("hashlib").sha256(summary.read_bytes()).hexdigest(),
        "n": n,
        "q_hat_normalization": rational(sum(h_raw)),
        "B_over_A_hat": rational(b),
        "E_qhat_g": rational(eg),
        "delta_BA_minus_Eg": rational(delta),
        "g_oscillation": rational(osc),
        "lower_bound_max_abs_BA_error_or_TV": rational(lower),
        "threshold_11_over_50": rational(threshold),
        "strictly_exceeds_11_over_50": lower > threshold,
        "inequality": "max(|B/A_hat-B/A_true|, TV(q_hat,q_true)) >= delta/(1+osc(g)) = 45*delta/67",
    }
    output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": "passed", "lower_bound": str(lower),
                      "exceeds_11_over_50": lower > threshold}))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("summary", type=Path); p.add_argument("output", type=Path)
    a = p.parse_args(); main(a.summary, a.output)
