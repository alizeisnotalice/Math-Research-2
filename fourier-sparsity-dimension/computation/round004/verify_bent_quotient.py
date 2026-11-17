#!/usr/bin/env python3
"""Independent exhaustive inverse-spectrum test for the 16-point bent sign class."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys


def transform(a):
    a = a[:]
    for step in (1, 2, 4, 8):
        for start in range(0, 16, 2 * step):
            for j in range(start, start + step):
                x, y = a[j], a[j + step]
                a[j], a[j + step] = x + y, x - y
    return a


def reduce_mod(x, basis):
    for b in basis:
        x = min(x, x ^ b)
    return x


parser = argparse.ArgumentParser()
parser.add_argument("result", type=Path)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
record = json.loads(args.result.read_text())['by_s']['16']
assert record['n'] == 4 and record['s'] == 16
base_signs = [1 if c > 0 else -1 for c in record['integercoeffs']]
# Echelon generators of constants and all four coordinate evaluations.
basis = [65535, 21845, 13107, 3855, 255]
valid = set()
for mask in range(1 << 16):
    coefficients = [base_signs[i] * (-1 if (mask >> i) & 1 else 1) for i in range(16)]
    if all(abs(v) == 4 for v in transform(coefficients)):
        valid.add(mask)
assert len(valid) == 896
quotient = {reduce_mod(x, basis) for x in valid}
assert len(quotient) == 28 and 0 in quotient
nonzero = sorted(quotient - {0})
pairs = [(a, b) for a, b in itertools.combinations(nonzero, 2) if (a ^ b) in quotient]
assert pairs
independent_triples = 0
for a, b, c in itertools.combinations(nonzero, 3):
    space = {0, a, b, c, a ^ b, a ^ c, b ^ c, a ^ b ^ c}
    if len(space) == 8:
        independent_triples += 1
        assert not space <= quotient
assert record['r'] == 2 and record['dimV'] == 7
# Check the saved graph's actual linear and affine ranks by complete value vectors.
coords = [sum(((i >> j) & 1) << i for i in range(16)) for j in range(4)]
def gf2_rank(rows):
    pivots = {}
    for x in rows:
        while x:
            p = x.bit_length() - 1
            if p not in pivots:
                pivots[p] = x
                break
            x ^= pivots[p]
    return len(pivots)
raw_rank = gf2_rank(coords + record['lifted_signcode_basis'])
lift_rank = gf2_rank(coords + record['lifted_signcode_basis'] + [65535])
assert raw_rank == 6 and lift_rank == 7
assert record['lifted_signcode_basis'][0] in valid
# A one-point sign error is illegal; this exercises the Booleanity rejection path.
damaged = record['lifted_signcode_basis'][0] ^ 1
assert damaged not in valid
result = {"status": "passed", "python": sys.version, "arithmetic": "integer",
          "input_sha256": hashlib.sha256(args.result.read_bytes()).hexdigest(),
          "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          "all_sign_patterns": 65536, "valid_sign_patterns": 896,
          "quotient_members": 28, "quotient_max_dimension": 2,
          "independent_triples_rejected": independent_triples,
          "raw_graph_linear_rank": raw_rank, "graph_affine_dimension": lift_rank - 1,
          "linear_rank_after_independent_character": lift_rank,
          "single_sign_corruption_rejected": True}
args.out.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result))
