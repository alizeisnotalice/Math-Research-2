#!/usr/bin/env python3
"""Exact residue-pair audit for the repaired L=10 middle-branch proof."""
import json

L = 10
A = (0, 1, 2, 5, 8)
B = (0, 2, 3)
rows = {r: [] for r in range(L)}
for i in A:
    for j in B:
        carry, r = divmod(i+j, L)
        factor = "v*" if carry else ""
        rows[r].append(f"{factor}A{i}B{j}")
expected = {
    0: ["A0B0", "v*A8B2"],
    1: ["A1B0", "v*A8B3"],
    2: ["A0B2", "A2B0"],
    3: ["A0B3", "A1B2"],
    4: ["A1B3", "A2B2"],
    5: ["A2B3", "A5B0"],
    6: [],
    7: ["A5B2"],
    8: ["A5B3", "A8B0"],
    9: [],
}
assert rows == expected, (rows, expected)
print(json.dumps({"status":"verified", "rows":rows}, indent=2))
