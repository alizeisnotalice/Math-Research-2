#!/usr/bin/env python3
"""Independent CPU validation of host-read terminal roots and GPU histogram."""
import json
from pathlib import Path

ROOT = Path(__file__).parent
P = 17


def replicas(path):
    return [json.loads(line) for line in path.read_text().splitlines()
            if json.loads(line).get("kind") == "replicate"]


def check_pair(host, base, label):
    assert len(host) == len(base), label
    summary = []
    for h, b in zip(host, base):
        for key in ("tau", "replicate", "seed_hex", "Z", "U", "logZ", "dead"):
            assert h[key] == b[key], (label, h["tau"], key, h[key], b[key])
        n, tau, N = h["_n"], h["tau"], h["_n"] ** 2
        roots = h["root_pairs"]
        allowed = h["terminal_allowed"]
        assert len(roots) == P and len(allowed) == P
        if h["dead"]:
            assert h["lastT"] == 0 and allowed == [0] * P
            assert h["hist_raw"] == [0] * (2 * n - 1)
            state = "dead"
        else:
            # Immediately before the last insertion the connected set has N-1
            # vertices, hence exactly one legal terminal vertex per particle.
            assert h["lastT"] == P and allowed == [1] * P
            if tau <= N:
                expected = [0] * (2 * n - 1)
                for r1, r2 in roots:
                    assert 0 <= r1 < N and 0 <= r2 < N
                    distance = abs(r1 % n - r2 % n) + abs(r1 // n - r2 // n)
                    expected[distance] += 1
                assert h["hist_raw"] == expected, (label, tau, h["hist_raw"], expected)
                assert sum(h["hist_raw"]) == h["lastT"] == P
                state = "alive-hist"
            else:
                assert h["hist_raw"] is None
                # no second root is selected in the one-peak convention
                assert all(r2 == 0 for _, r2 in roots)
                state = "alive-onepeak"
        summary.append({"n": n, "tau": tau, "state": state,
                        "lastT": h["lastT"], "hist_sum": None if h["hist_raw"] is None else sum(h["hist_raw"])})
    return summary


def add_n(rows, n):
    for row in rows:
        row["_n"] = n
    return rows


all_summary = []
all_summary += check_pair(add_n(replicas(ROOT / "n3-host-audit.jsonl"), 3),
                          replicas(ROOT / "n3-base.jsonl"), "n3-all-taus")
for n, taus in ((16, (2, 128, 255, 256, 257)),
                (24, (2, 288, 575, 576, 577)),
                (32, (2, 512, 1023, 1024, 1025))):
    for tau in taus:
        all_summary += check_pair(add_n(replicas(ROOT / f"host-n{n}-tau{tau}.jsonl"), n),
                                  replicas(ROOT / f"base-n{n}-tau{tau}.jsonl"),
                                  f"n{n}-tau{tau}")

certificate = {
    "status": "passed",
    "method": "independent Python Manhattan bins from host-read terminal State roots; terminal allowed and frozen-base terminal-field comparisons",
    "particles": P,
    "cases": all_summary,
    "fk_normalization": {
        "theory": "For an alive run, H_b is the unweighted terminal particle count and t_N=lastT=P because all P preterminal connected N-1-vertex states have one final legal insertion. The bin FK estimator is Z*(H_b/P)=Z*(H_b/lastT).",
        "dead_case": "If a preceding rank has zero total weight, lastT=0, the run is dead, and its all-zero histogram is not divided by lastT.",
    },
}
(ROOT / "histogram-independent-certificate.json").write_text(json.dumps(certificate, indent=2) + "\n")
print(json.dumps(certificate, indent=2))
