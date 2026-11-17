#!/usr/bin/env python3
"""Audit derived first-merge output. Finite diagnostics, not a mixing certificate."""
import argparse
import hashlib
import json
import math
from pathlib import Path

DERIVED = "3a9f3c0a8c10409346c25f16f2d82463f928bfaacf7307b4a88a770ec8f69ec6"
MIN_SUB = float.fromhex("0x0.0000000000001p-1022")

def require(ok, message):
    if not ok:
        raise ValueError(message)

def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)

def close(x, y):
    return math.isclose(x, y, rel_tol=2e-10, abs_tol=2e-11)

def unique_records(records, kind, count):
    rows = [r for r in records if r.get("kind") == kind]
    require(len(rows) == count, kind + " record count")
    result = {}
    for row in rows:
        i = row.get("replicate")
        require(type(i) is int and 0 <= i < count and i not in result, kind + " duplicate/invalid ID")
        result[i] = row
    return result

def audit(raw, source, exact_dir=None):
    require(hashlib.sha256(source.read_bytes()).hexdigest() == DERIVED, "source mismatch")
    records = [json.loads(s) for s in raw.read_text().splitlines() if s.strip()]
    require(len(records) >= 4, "incomplete output")
    meta, end = records[0], records[-1]
    require(meta.get("kind") == "metadata" and end.get("kind") == "run_end", "missing metadata/run_end")
    require(meta.get("implementation") == "r026-tau2-firstmerge-rb-terminal", "implementation mismatch")
    n, P, R, L = [meta.get(k) for k in ("n", "particles", "replicates", "move_steps")]
    require(all(type(v) is int for v in (n,P,R,L)), "noninteger metadata")
    require(3 <= n <= 32 and P > 0 and R > 0 and L >= 0, "metadata range")
    require(meta.get("tau_lo") == meta.get("tau_hi") == 2 and meta.get("method") in ("strat", "multi"), "target/method")
    require(meta.get("rb_payload_bytes_per_particle") == 8 and meta.get("rb_gpu_working_bytes_per_particle") == 256, "payload layout")
    require(len(records) == 2 * R + 2, "extra/missing records")
    require(all(r.get("kind") in ("replicate", "firstmerge_rb_calibration") for r in records[1:-1]), "unexpected record kind")
    require(finite(end.get("wall_seconds")) and end["wall_seconds"] >= 0, "end timing")
    domain = int(meta["domain"], 16)
    require(0 <= domain < 2**64, "domain range")
    reps = unique_records(records, "replicate", R)
    rb = unique_records(records, "firstmerge_rb_calibration", R)
    group_count = (P + 255) // 256
    for i in range(R):
        x, y = reps[i], rb[i]
        require(x.get("tau") == y.get("tau") == 2 and x.get("dead") is False, "target/dead particle")
        seed = (0xA022000000000000 + domain + n*1000000 + 20000 + i) % 2**64
        require(int(x["seed_hex"],16) == seed, "seed mismatch")
        require(finite(x.get("logZ")) and x["logZ"] <= 1e-10 and finite(x.get("Z")) and x["Z"] >= 0, "invalid normalizer")
        expected_z = math.exp(x["logZ"])
        require(math.isclose(x["Z"], expected_z, rel_tol=5e-12, abs_tol=MIN_SUB), "Z/logZ inconsistent")
        for key in ("mean_m", "mean_m_d2", "mean_m_F0", "mean_m_G0"):
            require(finite(y.get(key)), "nonfinite RB field")
        require(1-1e-10 <= y["mean_m"] <= (n*n-1)/2+1e-10 and y["mean_m_d2"] > 0 and y["mean_m_F0"] >= 0, "RB range")
        h = y["mean_m_hist"]
        require(len(h) == 2*n-1 and all(finite(v) and v >= 0 for v in h) and h[0] == h[1] == 0, "RB histogram")
        require(close(math.fsum(h), y["mean_m"]), "RB histogram mass")
        require(close(y["mean_m_F0"]-2*y["mean_m_G0"], 4*y["mean_m_d2"]), "RB geometry identity")
        groups = [x[k] for k in ("tau2_lift_L_groups", "tau2_lift_L_distance2_groups", "tau2_lift_L_F0_groups", "tau2_lift_L_G0_groups")]
        require(all(len(a) == group_count and all(type(v) is int for v in a) for a in groups), "integer group layout")
        for a,b,c,d in zip(*groups):
            require(a > 0 and b > 0 and c >= 0 and c-2*d == 4*b, "integer group identity")
        hist = x["tau2_distance_L_hist_groups"]
        require(len(hist) == group_count, "L histogram groups")
        for j, h in enumerate(hist):
            require(len(h) == 2*n-1 and all(type(v) is int and v >= 0 for v in h) and h[0] == h[1] == 0 and sum(h) == groups[0][j], "L histogram identity")
        keys = ("tau2_lift_scaled_mean_L", "tau2_lift_scaled_mean_L_distance2", "tau2_lift_scaled_mean_L_F0", "tau2_lift_scaled_mean_L_G0")
        for key, group in zip(keys, groups):
            require(finite(x.get(key)) and close(x[key], sum(group)/P), "L mean mismatch")
    mx = max(x["logZ"] for x in reps.values())
    w = [math.exp(reps[i]["logZ"]-mx) for i in range(R)]
    s = math.fsum(w)
    av = lambda key: math.fsum(w[i]*rb[i][key] for i in range(R))/s
    mh = av("mean_m")
    logb = mx + math.log(s/R)
    lifted = [w[i]*rb[i]["mean_m"] for i in range(R)]
    sl = math.fsum(lifted)
    lmean = math.fsum(w[i]*reps[i]["tau2_lift_scaled_mean_L"] for i in range(R))/s
    ld = math.fsum(w[i]*reps[i]["tau2_lift_scaled_mean_L_distance2"] for i in range(R))/s
    out = {"audit_passed": True, "audit_scope": "input and implementation invariants only; no SMC accuracy or mixing certification", "derived_source_sha256": DERIVED, "raw_sha256": hashlib.sha256(raw.read_bytes()).hexdigest(), "n": n, "P": P, "R": R, "L": L, "domain": meta["domain"], "method": meta["method"], "ESS": s*s/math.fsum(v*v for v in w), "maxshare": max(w)/s, "m_lift_ESS": sl*sl/math.fsum(v*v for v in lifted), "m_lift_maxshare": max(lifted)/sl, "log_B_over_Nfactorial": logb, "log_mA2_over_Nfactorial": logb+math.log(mh), "B_over_Nfactorial_display": math.exp(logb), "mA2_over_Nfactorial_display": math.exp(logb+math.log(mh)), "scaled_mean_m": mh, "scaled_mean_mD": av("mean_m_d2"), "scaled_mean_mF": av("mean_m_F0"), "scaled_mean_mG": av("mean_m_G0"), "A2_over_B_mRB": mh, "Q_mRB": av("mean_m_d2")/(n*n*mh), "Q_L_legacy": ld/(n*n*lmean), "distance_mRB": [math.fsum(w[i]*rb[i]["mean_m_hist"][d] for i in range(R))/(s*mh) for d in range(2*n-1)]}
    if exact_dir is not None and (exact_dir/f"n{n}.json").exists():
        path = exact_dir/f"n{n}.json"
        exact = json.loads(path.read_text())
        require(exact.get("n") == n, "exact size mismatch")
        A, B = int(exact["A2"]), int(exact["Btop2"])
        require(0 < B <= A <= math.factorial(n*n), "exact counts range")
        out.update(exact_A2_over_B=A/B, exact_source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    else:
        out["exact_comparison"] = "unavailable: no matching exact input"
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument("raw", type=Path)
    p.add_argument("source", type=Path)
    p.add_argument("--exact-dir", type=Path)
    args = p.parse_args()
    print(json.dumps(audit(args.raw,args.source,args.exact_dir), indent=2, allow_nan=False))

if __name__ == "__main__":
    main()
