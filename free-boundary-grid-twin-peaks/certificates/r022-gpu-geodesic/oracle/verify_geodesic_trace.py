#!/usr/bin/env python3
"""Recompute every particle's terminal prefix flags and 256-group sums."""
import json
import math
import pathlib
import sys

if len(sys.argv) not in (2, 3):
    raise SystemExit("usage: verify_geodesic_trace.py trace.json [run.jsonl]")

trace_path = pathlib.Path(sys.argv[1])
trace = json.loads(trace_path.read_text())
n, P = trace["n"], trace["P"]
m = 2 * n - 1
terminal = next(row for row in trace["rows"] if row.get("terminal"))
prefixes = terminal["prefix"]
assert len(prefixes) == P
assert all(len(word) >= m for word in prefixes), (n, m, min(map(len, prefixes)))

def md(a, b):
    return abs(a % n - b % n) + abs(a // n - b // n)

def indicator(word, d):
    return (all(md(word[k - 1], word[k]) == 1 for k in range(1, d))
            and md(word[0], word[d - 1]) == d - 1)

# This catches the tempting but incorrect endpoint-only implementation:
# for n>=3, 0 -> 4 -> 2 has endpoints at distance two but has two jumps.
if n >= 3:
    assert indicator([0, 1, 2], 3)
    assert not indicator([0, n + 1, 2], 3)

expected = [[0] * m for _ in range((P + 255) // 256)]
for particle, word in enumerate(prefixes):
    assert len(set(word)) == len(word)
    assert all(0 <= x < n * n for x in word)
    for d in range(1, m + 1):
        expected[particle // 256][d - 1] += int(indicator(word, d))

actual = [[int(x) for x in group]
          for group in terminal["geodesic_prefix_groups"]]
assert actual == expected, {"expected": expected, "actual": actual}

result = {"passed": True, "n": n, "P": P, "m": m,
          "groups": expected,
          "raw": [sum(group[k] for group in expected) for k in range(m)]}
if len(sys.argv) == 3:
    records = [json.loads(line) for line in pathlib.Path(sys.argv[2]).read_text().splitlines()]
    rep = next(x for x in records if x["kind"] == "replicate")
    assert rep["geodesic_prefix_group_counts"] == expected
    denom = int(terminal["cdf"][-1])
    want = [x / denom for x in result["raw"]]
    got = rep["scaled_geodesic_prefix_p"]
    assert len(got) == m
    assert all(math.isclose(a, b, rel_tol=0.0, abs_tol=0.0)
               for a, b in zip(got, want)), {"got": got, "want": want}
    result["terminal_total"] = denom
    result["scaled"] = want
print(json.dumps(result, indent=2))
