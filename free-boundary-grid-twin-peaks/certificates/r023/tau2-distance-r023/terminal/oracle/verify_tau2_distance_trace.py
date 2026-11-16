#!/usr/bin/env python3
"""Independently reconstruct r023 L-weighted root-distance fields from trace."""
import json
import pathlib
import sys

if len(sys.argv) not in (2, 3):
    raise SystemExit("usage: verify_tau2_distance_trace.py trace.json [run.jsonl]")

trace = json.loads(pathlib.Path(sys.argv[1]).read_text())
n, P, N = trace["n"], trace["P"], trace["n"] ** 2
m, G = 2 * n - 1, (P + 255) // 256
terminal = next(row for row in trace["rows"] if row.get("terminal"))
assert terminal["tau2_distance_bad"] == 0
assert terminal["tau2_lift_bad"] == 0


def md(a, b):
    return abs(a % n - b % n) + abs(a // n - b // n)


hist = [[0] * m for _ in range(G)]
l2 = [0] * G
for particle, word in enumerate(terminal["prefix"]):
    assert len(word) == N - 1 and len(set(word)) == N - 1
    state = terminal["state"][particle]
    vx, vy = state[2], state[3]
    assert word[:2] == [vx, vy] and md(vx, vy) != 1
    best = min(rank + 1 for rank, vertex in enumerate(word) if md(vy, vertex) == 1)
    assert 2 < best < N
    L, d = best - 2, md(vx, vy)
    assert 1 <= L <= N - 3 and 2 <= d < m
    hist[particle // 256][d] += L
    l2[particle // 256] += L * L

assert all(group[0] == group[1] == 0 for group in hist)
lift_l = [int(x) for x in terminal["tau2_lift_L_groups"]]
assert [sum(group) for group in hist] == lift_l
assert hist == [[int(x) for x in group] for group in terminal["tau2_distance_L_hist_groups"]]
assert l2 == [int(x) for x in terminal["tau2_distance_L2_groups"]]
result = {"passed": True, "n": n, "P": P, "hist_groups": hist, "L2_groups": l2,
          "raw_hist": [sum(group[k] for group in hist) for k in range(m)],
          "raw_L2": sum(l2)}
if len(sys.argv) == 3:
    records = [json.loads(line) for line in pathlib.Path(sys.argv[2]).read_text().splitlines()]
    rep = next(row for row in records if row["kind"] == "replicate")
    assert rep["tau2_distance_L_hist_groups"] == hist
    assert rep["tau2_distance_L2_groups"] == l2
    expected_hist = [value / P for value in result["raw_hist"]]
    assert rep["tau2_distance_scaled_mean_L_hist"] == expected_hist
    assert rep["tau2_distance_scaled_mean_L2"] == result["raw_L2"] / P
print(json.dumps(result, indent=2))
