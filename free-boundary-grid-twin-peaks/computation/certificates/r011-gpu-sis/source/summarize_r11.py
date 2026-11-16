#!/usr/bin/env python3
"""Aggregate r11 xoshiro256** JSONL repetitions.

The per-repetition delta diagnostic uses only within-tau sampling variance:
 Var(sum_tau mean(W(f-NQ))) = sum_tau [sum x^2-(sum x)^2/M]/[M(M-1)].
It deliberately excludes between-tau mean differences and is not a certificate.
"""
import glob, json, math, statistics
from pathlib import Path


def one(path):
    rows = [json.loads(s) for s in Path(path).read_text().splitlines()]
    summary = next(x for x in rows if x["kind"] == "summary")
    n, q, L = summary["n"], summary["Q_euclidean_hat"], summary["log_A2_hat"]
    var = 0.0
    for x in rows:
        if x["kind"] != "stratum" or x["alive"] == 0:
            continue
        m, c = x["samples"], n * n * q
        e1 = math.exp(x["log_EW"] - L)
        ef = math.exp(x["log_EWf"] - L)
        e2 = math.exp(x["log_EW2"] - 2 * L)
        e2f = math.exp(x["log_EW2f"] - 2 * L)
        e2f2 = math.exp(x["log_EW2f2"] - 2 * L)
        sx = m * (ef - c * e1)
        sx2 = m * (e2f2 - 2 * c * e2f + c * c * e2)
        var += (sx2 - sx * sx / m) / (m * (m - 1))
    out = dict(path=path, q=q, delta_se=math.sqrt(max(var, 0)) / (n * n),
               gpu_seconds=summary["gpu_seconds"], wall_seconds=summary["wall_seconds"])
    return out


if __name__ == "__main__":
    reps = [one(p) for p in sorted(glob.glob("runs/n8-rep*-r11/results.jsonl"))]
    qs = [x["q"] for x in reps]
    result = {"repetitions": reps, "mean_q": statistics.mean(qs),
              "sample_sd_q": statistics.stdev(qs),
              "sem_across_reps": statistics.stdev(qs) / math.sqrt(len(qs)),
              "warning": "floating-point deterministic-RNG numerical diagnostic, not a confidence certificate"}
    Path("summary-r11-n8.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
