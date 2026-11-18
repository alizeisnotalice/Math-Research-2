"""Verify target-sieve factor degree claims from Frobenius polynomial gcds."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from math import lcm
from pathlib import Path

from sympy import ZZ
from sympy.polys.galoistools import gf_gcd, gf_pow_mod, gf_sub


H = (1, -4, 6, -7, 4, -2, 1)


def certify(item):
    coefficients, prime, reported = item
    expected = Counter(degree for degree, _ in reported)
    polynomial = list(coefficients)
    x = [1, 0]
    power = x
    found = {}
    # Compare the full gcd-degrees, not just the lcm from factor_list.
    for degree in range(1, len(polynomial)):
        power = gf_pow_mod(power, prime, polynomial, prime, ZZ)
        difference = gf_sub(power, x, prime, ZZ)
        gcd_degree = len(gf_gcd(polynomial, difference, prime, ZZ))-1
        lower = sum(d*count for d, count in found.items() if degree % d == 0)
        remaining = gcd_degree-lower
        assert remaining >= 0 and remaining % degree == 0
        count = remaining//degree
        assert count == expected[degree], (prime, degree, count, expected[degree])
        if count:
            found[degree] = count
    return lcm(*found)


def main():
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    raw = args.packet.read_bytes()
    packet = json.loads(raw)
    from sympy import Poly, Symbol
    x = Symbol("x")
    jobs = []
    required = {int(p): int(v) for p,v in packet["h_reduction_periods"].items()}
    expected = []
    for p, period in required.items():
        factors = Poly.from_list(H,x,modulus=p).factor_list()[1]
        jobs.append((H,p,[(f.degree(),e) for f,e in factors]))
        expected.append(period)
    seen = set()
    for row in packet["rejected_draws"]:
        assert row["kind"] == "exact_residue_period"
        prime = row["prime"]
        degrees = row["degrees"]
        period = lcm(*(degree for degree,_ in degrees))
        assert period % required[prime]
        assert row["required"] == required[prime]
        assert row["draw"] not in seen
        seen.add(row["draw"])
        jobs.append((tuple(reversed(row["C"])), prime, degrees))
        expected.append(period)
    for task in packet["tasks"]:
        assert task["draw"] not in seen
        seen.add(task["draw"])
        assert task["C"][:12] == [1]*7+[0]+[1]*4
        assert set(map(int,task["periods"])) == set(required)
        for p in required:
            row = task["periods"][str(p)]
            assert row["period"] % required[p] == 0
            jobs.append((tuple(reversed(task["C"])),p,row["degrees"]))
            expected.append(row["period"])
    assert seen == set(range(packet["attempted_draws"]))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        assert list(pool.map(certify, jobs)) == expected
    report={"status":"verified","packet_sha256":sha256(raw).hexdigest(),"period_tests":len(jobs),
            "exact_rejected":len(packet["rejected_draws"]),"numerical_targets":len(packet["tasks"])}
    with args.out.open("x") as output:
        json.dump(report,output,indent=2)
        output.write("\n")
    print(json.dumps(report))


if __name__=="__main__":
    main()
