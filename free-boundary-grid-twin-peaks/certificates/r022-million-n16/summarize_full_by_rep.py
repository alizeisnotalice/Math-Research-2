#!/usr/bin/env python3
"""Aggregate a multi-tau SMC JSONL run by independent replicate.

Each tau uses its own particle population but shares a replicate index.  The
unbiased full normalizer estimator for repeat r is therefore
    sum_tau Z[tau, r],
not the pool of all tau/repeat records.  This script reads finite logZ values
only (a dead particle system contributes exactly zero), uses Decimal-80 for
all exponentiation/sums, and rejects incomplete tau x replicate rectangles.
The output is a numerical summary, not a confidence certificate.
"""
import json
import sys
from collections import defaultdict
from decimal import Decimal, getcontext

getcontext().prec = 80


def D(x):
    return Decimal(str(x))


def stable_weights(logs):
    """Return (sum exp(logs), normalized weights) without float arithmetic."""
    if not logs:
        return Decimal(0), []
    m = max(logs)
    w = [(x - m).exp() for x in logs]
    s = sum(w)
    return m.exp() * s, [x / s for x in w]


def scalar_summary(values):
    R = len(values)
    mean = sum(values) / Decimal(R)
    var = (sum((x - mean) ** 2 for x in values) / Decimal(R - 1)
           if R > 1 else Decimal(0))
    return mean, var, (var / Decimal(R)).sqrt()


def text(x):
    return str(x) if x is not None else None


def summarize_tau(records):
    logs = [D(x["logZ"]) for x in records if not x["dead"]]
    z, shares = stable_weights(logs)
    zvals = [(D(x["logZ"]).exp() if not x["dead"] else Decimal(0))
             for x in records]
    mean, var, se = scalar_summary(zvals)
    live = [x for x in records if not x["dead"]]
    numerator = sum(D(x["logZ"]).exp() * D(x["U_over_Z"])
                    for x in live if x.get("U_over_Z") is not None)
    denominator = sum(D(x["logZ"]).exp()
                      for x in live if x.get("U_over_Z") is not None)
    qs = [D(x["U_over_Z"]) for x in live if x.get("U_over_Z") is not None]
    return {
        "tau": records[0]["tau"], "R": len(records),
        "dead": len(records) - len(live),
        "mean_Z_decimal": text(mean),
        "sample_variance_Z_decimal": text(var),
        "descriptive_SE_Z_decimal": text(se),
        "repeat_logZ_ESS": text(Decimal(1) / sum(w * w for w in shares)) if shares else "0",
        "repeat_max_share": text(max(shares)) if shares else "0",
        "logZ_min": text(min(logs)) if logs else None,
        "logZ_max": text(max(logs)) if logs else None,
        "Z_weighted_U_over_Z": text(numerator / denominator) if denominator else None,
        "unweighted_mean_U_over_Z": text(sum(qs) / Decimal(len(qs))) if qs else None,
    }


def main(path):
    rows = [json.loads(line) for line in open(path, encoding="utf-8") if line.strip()]
    meta = next((x for x in rows if x.get("kind") == "metadata"), None)
    if meta is None:
        raise SystemExit("missing metadata")
    R = int(meta["replicates"])
    lo, hi = int(meta["tau_lo"]), int(meta["tau_hi"])
    records = [x for x in rows if x.get("kind") == "replicate"]
    by_tau = defaultdict(list)
    by_rep = defaultdict(list)
    for x in records:
        tau, rep = int(x["tau"]), int(x["replicate"])
        if not lo <= tau <= hi or not 0 <= rep < R:
            raise SystemExit(f"out-of-range record tau={tau}, replicate={rep}")
        by_tau[tau].append(x)
        by_rep[rep].append(x)
    expected_taus = set(range(lo, hi + 1))
    if set(by_tau) != expected_taus or set(by_rep) != set(range(R)):
        raise SystemExit("incomplete tau or replicate set")
    for tau in expected_taus:
        if len(by_tau[tau]) != R or {x["replicate"] for x in by_tau[tau]} != set(range(R)):
            raise SystemExit(f"incomplete or duplicate tau={tau}")
    for rep in range(R):
        if len(by_rep[rep]) != len(expected_taus) or {x["tau"] for x in by_rep[rep]} != expected_taus:
            raise SystemExit(f"incomplete or duplicate replicate={rep}")

    tau_summary = [summarize_tau(sorted(by_tau[t], key=lambda x: x["replicate"]))
                   for t in range(lo, hi + 1)]
    total_z, total_u, total_logz, all_dead = [], [], [], 0
    for rep in range(R):
        rr = sorted(by_rep[rep], key=lambda x: x["tau"])
        live = [x for x in rr if not x["dead"]]
        logs = [D(x["logZ"]) for x in live]
        z, _ = stable_weights(logs)
        if not live:
            all_dead += 1
            total_z.append(Decimal(0)); total_u.append(Decimal(0)); continue
        total_z.append(z)
        total_logz.append(z.ln())
        total_u.append(sum(D(x["logZ"]).exp() * D(x["U_over_Z"])
                           for x in live if x.get("U_over_Z") is not None))
    mean, var, se = scalar_summary(total_z)
    mean_u = sum(total_u) / Decimal(R)
    _, shares = stable_weights(total_logz)
    output = {
        "kind": "multi_tau_by_replicate_summary",
        "input": path,
        "metadata": meta,
        "records": len(records),
        "expected_records": (hi - lo + 1) * R,
        "full_replicates_all_dead": all_dead,
        "full_mean_sum_tau_Z_decimal": text(mean),
        "full_sample_variance_sum_tau_Z_decimal": text(var),
        "full_descriptive_SE_sum_tau_Z_decimal": text(se),
        "full_log_mean_sum_tau_Z": text(mean.ln()) if mean else None,
        "full_repeat_logZ_ESS": text(Decimal(1) / sum(w * w for w in shares)) if shares else "0",
        "full_repeat_max_share": text(max(shares)) if shares else "0",
        "full_Z_weighted_U_over_Z": text(mean_u / mean) if mean else None,
        "per_tau": tau_summary,
    }
    print(json.dumps(output, sort_keys=True))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} RUN.jsonl")
    main(sys.argv[1])
