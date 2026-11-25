"""Dual-basis lower bound on the length of every nonzero integral relation."""

import json
from fractions import Fraction
from flint import fmpz_mat, fmpq_mat
from lattice import block


def audit(n):
    data = block(n)
    V = fmpz_mat([v["v"] for v in data["vectors"]])
    G = V * V.transpose()
    inverse = fmpq_mat(G).inv()
    worst = max(
        Fraction(int(inverse[j, j].p), int(inverse[j, j].q)) for j in range(V.nrows())
    )
    lower = 1 / worst
    upper = min(v["length"] for v in data["vectors"])
    return dict(
        N=n,
        nullity=data["nullity"],
        image_rank=data["image_rank"],
        lower_squared_floor=lower.numerator // lower.denominator,
        upper_l1=upper,
        excluded=lower > 2**n,
    )


if __name__ == "__main__":
    print(json.dumps([audit(n) for n in (32, 40, 48, 64)], indent=2))
