#!/usr/bin/env python3
"""Independent integer verification of every recorded round004 sign-code witness."""
import argparse
import hashlib
import json
from pathlib import Path

def fwht(values):
    values = values[:]
    step = 1
    while step < len(values):
        for start in range(0, len(values), 2 * step):
            for j in range(start, start + step):
                a, b = values[j], values[j + step]
                values[j], values[j + step] = a + b, a - b
        step *= 2
    return values

def rank(vectors):
    basis = {}
    for v in vectors:
        x = v
        while x:
            bit = x.bit_length() - 1
            if bit in basis:
                x ^= basis[bit]
            else:
                basis[bit] = x
                break
    return len(basis)

def verify(row):
    n, support, coefficients, lifts = row["n"], row["support"], row["integercoeffs"], row["lifted_signcode_basis"]
    if len(support) != row["s"] or len(coefficients) != row["s"]:
        raise AssertionError("support length")
    if rank(row["C_basis"] + lifts) != row["dimV"] or len(lifts) != row["r"]:
        raise AssertionError("witness rank")
    scale = 1 << n
    for z in range(1 << len(lifts)):
        walsh = [0] * scale
        for i, frequency in enumerate(support):
            exponent = sum(((z >> j) & 1) * ((lifts[j] >> i) & 1) for j in range(len(lifts))) & 1
            walsh[frequency] = -coefficients[i] if exponent else coefficients[i]
        truth = fwht(walsh)
        if any(abs(value) != scale for value in truth):
            raise AssertionError((row["n"], row["s"], z, truth))
    return {"n": n, "s": row["s"], "dimV": row["dimV"], "r": row["r"], "slices": 1 << len(lifts)}

parser = argparse.ArgumentParser()
parser.add_argument("result", type=Path)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
raw = args.result.read_bytes()
data = json.loads(raw)
rows = [row for per_n in data["by_n"].values() for row in per_n["rows"]]
checked = [verify(row) for row in rows]
args.out.write_text(json.dumps({"status": "passed", "input_sha256": hashlib.sha256(raw).hexdigest(),
                                "rows_verified": len(checked), "total_boolean_slices": sum(x["slices"] for x in checked),
                                "max_dimV": max(x["dimV"] for x in checked)}, indent=2) + "\n")
