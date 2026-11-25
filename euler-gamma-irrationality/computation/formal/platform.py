"""Certificate of the neighboring fractional-part comparison on every plateau."""

import json
import sys
from collections import Counter
from math import lcm
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from intervals import certify
from sondow import case


def run(end):
    cases = {}
    results = []
    for n in range(1, end + 1):
        d = lcm(*range(1, 2 * n + 1))
        if lcm(d, 2 * n + 1, 2 * n + 2) != d:
            continue
        for i in (n, n + 1):
            if i not in cases:
                cases[i] = case(i)
        left, right = cases[n], cases[n + 1]
        coeff = Counter({p: 16 * v for p, v in right["prime_exponents"].items()})
        coeff.subtract(left["prime_exponents"])
        j, k = left["log_floor"]["decision"], right["log_floor"]["decision"]
        precision = sum(map(abs, coeff.values())).bit_length() + 48
        bound = certify(0, j - 16 * k, coeff, start=precision, limit=precision + 256)
        assert bound["decision"] is not None
        results.append(dict(n=n, d=d, floors=[j, k], difference=bound))
    return results


if __name__ == "__main__":
    print(json.dumps(run(128), indent=2))
