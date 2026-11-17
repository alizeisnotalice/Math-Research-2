#!/usr/bin/env python3
"""Row-transfer resource estimate for the selected-prefix IE kernel (no C(I))."""
from collections import defaultdict
from math import factorial
import json
from pathlib import Path


def independent_histogram(rows: int, cols: int) -> dict[int, int]:
    valid = [m for m in range(1 << cols) if not (m & (m << 1))]
    cur = {(0, 0): 1}
    for _ in range(rows):
        nxt = defaultdict(int)
        for (previous, size), multiplicity in cur.items():
            for mask in valid:
                if not (previous & mask):
                    nxt[mask, size + mask.bit_count()] += multiplicity
        cur = nxt
    out = defaultdict(int)
    for (_, size), multiplicity in cur.items():
        out[size] += multiplicity
    return dict(out)


def record(prefix: str, rows: int) -> dict:
    hist = independent_histogram(rows, 7)
    count = sum(hist.values())
    return {
        "prefix": prefix,
        "Uset_geometry": f"{rows}x7 rectangle",
        "independent_set_histogram": {str(k): hist[k] for k in sorted(hist)},
        "independent_sets": count,
        "max_independent_size": max(hist),
        "C_I_upper_bit_length": factorial(42).bit_length(),
        "signed_accumulator_coarse_upper_bit_length": (count * factorial(42)).bit_length(),
    }


out = {
    "n": 7,
    "runs": [record("top row S=0x7f", 5), record("first two rows S=0x3fff", 4)],
    "scope": "row-transfer cardinality estimate only; no n=7 C(I) or completion count is evaluated",
}
Path("n7-resource-estimate.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
