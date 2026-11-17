#!/usr/bin/env python3
"""Reconstruct four modular layer passes and compare frozen finite values."""
from __future__ import annotations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).parent
RUNS = ROOT / (sys.argv[1] if len(sys.argv) == 2 else "modular-runs-v2")
PRIMES = [9223372036854775783, 9223372036854775643, 9223372036854775549, 9223372036854775507]
EXPECTED = {
    3: 156,
    4: 4213440,
    5: 153450311313600,
    6: 20779014929894865426585600,
    7: 23386940526059745800104592446335589635840,
}

def crt(residues):
    x, modulus = 0, 1
    for r, p in zip(residues, PRIMES):
        t = ((r - x) % p) * pow(modulus % p, -1, p) % p
        x += modulus * t
        modulus *= p
    return x, modulus

out = []
for n, expected in EXPECTED.items():
    files = [RUNS / f"n{n}-p{i}.json" for i in range(4)]
    rows = [json.loads(p.read_text()) for p in files]
    assert [int(x["prime"]) for x in rows] == PRIMES
    residues = [int(x["residue"]) for x in rows]
    answer, modulus = crt(residues)
    assert all(r == expected % p for r, p in zip(residues, PRIMES))
    assert answer == expected and answer < modulus
    out.append({"n": n, "reconstructed_H_toprow": str(answer),
                "all_residues_match_exact": True, "crt_unique": True,
                "independent_sets": rows[0]["independent_sets"]})
print(json.dumps({"all_match": True, "rows": out,
                  "run_directory": RUNS.name, "crt_modulus": str(modulus), "crt_modulus_bit_length": modulus.bit_length()}, indent=2))
