#!/usr/bin/env python3
"""Summarize the one-peak probability estimator without confusing it with A1."""
import hashlib
import json
import math
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 90
HERE = Path(__file__).resolve().parent
RAW = HERE / "n12-p65536-r128-t145-l16-domain9300000000000000.jsonl"
TEL = HERE / "n12-p65536-r128-t145-l16.telemetry.jsonl"
N, P, R = 144, 65536, 128
M2MAX, M3MAX = 1020, 259080

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
records = [json.loads(x) for x in RAW.read_text().splitlines()]
rows = [x for x in records if x.get("kind") == "replicate"]
end = next(x for x in records if x.get("kind") == "run_end")
assert len(rows) == R and all(not x["dead"] and x["tau"] == N + 1 for x in rows)
assert all(x["K2"] is None for x in rows)
for x in rows:
    assert len(x["inverse_m2_groups"]) == (P + 255) // 256
    assert len(x["inverse_m3_groups"]) == (P + 255) // 256
    assert max(x["inverse_m2_groups"]) <= 256 * M2MAX
    assert max(x["inverse_m3_groups"]) <= 256 * M3MAX
    assert sum(x["inverse_m2_groups"]) <= P * M2MAX
    assert sum(x["inverse_m3_groups"]) <= P * M3MAX
sum_z = sum(Decimal(str(x["Z"])) for x in rows)
sum_zm2 = sum(Decimal(str(x["Zmean_m2"])) for x in rows)
sum_zm3 = sum(Decimal(str(x["Zmean_m3"])) for x in rows)
mean_z, mean_zm2, mean_zm3 = sum_z / R, sum_zm2 / R, sum_zm3 / R
fac = math.factorial(N)
facd = Decimal(fac)
tel = [json.loads(x) for x in TEL.read_text().splitlines()]
samples = [x for x in tel if x.get("kind") == "sample"]
rss = [x["rss_kib"] for x in samples if x.get("rss_kib") is not None]
agx = {k: [x.get(k) for x in samples if x.get(k) is not None]
       for k in ["Device Utilization %", "Renderer Utilization %", "Alloc system memory", "In use system memory"]}
out = {
 "status": "completed", "scope": "finite SMC numerical diagnostic; no exact n=12 count asserted",
 "raw_sha256": sha(RAW), "telemetry_sha256": sha(TEL),
 "run": {"n": 12, "N": N, "P": P, "R": R, "tau": N + 1, "L": 16, "domain": "9300000000000000",
         "run_wall_seconds": end["wall_seconds"], "sum_gpu_seconds": str(sum(Decimal(str(x["gpu_seconds"])) for x in rows))},
 "probability_normalization": {
   "meaning_of_Z": "A1/N!, the uniform-permutation probability of a one-peak word",
   "mean_Z": str(mean_z), "mean_Zmean_m2": str(mean_zm2), "mean_Zmean_m3": str(mean_zm3),
   "R2_estimate_A2d2_over_A1": str(sum_zm2 / sum_z), "R3_estimate_A2d3_over_A1": str(sum_zm3 / sum_z),
   "N_factorial_exact": str(fac),
   "A1_estimate_decimal": str(mean_z * facd), "A2d2_estimate_decimal": str(mean_zm2 * facd),
   "A2d3_estimate_decimal": str(mean_zm3 * facd),
   "log_A1_estimate": str(mean_z.ln() + facd.ln())
 },
 "integer_checks": {"per_word_m2_max": M2MAX, "per_word_m3_max": M3MAX,
                    "per_group_m2_max": 256*M2MAX, "per_group_m3_max": 256*M3MAX,
                    "per_replicate_m3_raw_max": P*M3MAX, "passed": True},
 "telemetry": {"scope": "target RSS; AGX device-wide and nonexclusive", "samples": len(samples),
               "rss_kib_min": min(rss) if rss else None, "rss_kib_max": max(rss) if rss else None,
               "agx_min_max": {k: ([min(v), max(v)] if v else None) for k,v in agx.items()},
               "exit": next(x for x in tel if x.get("kind") == "exit")}
}
(HERE / "n12-p65536-r128-t145-l16.summary.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
