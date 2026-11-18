"""Derive the seven-row elimination and check its exact residue collisions."""

import argparse
from collections import defaultdict
from itertools import product
import json
from pathlib import Path

import sympy as sp


def main():
    if not __debug__:
        raise RuntimeError("run with assertions enabled")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    t, u = sp.symbols("t u")
    h = t**6 - 4*t**5 + 6*t**4 - 7*t**3 + 4*t**2 - 2*t + 1

    def reduce(expression):
        return sp.rem(expression, h, t)

    m = 1 - t + t**2
    n = 1 - t*m
    p = 1 - (1-t)*n
    v = t**2*m
    a = (1, t, 0, n, p, 0, 0)
    b = (1, 1-t, m, 0, 0, v, 0)
    assert all(reduce(sum(a[i]*b[k-i] for i in range(k+1))) == 1 for k in range(7))

    # The eliminated coordinates are B0,B1,B2,A3,A4,B5 within
    # the seven-row block. Their response columns are shifts of a or b.
    columns = (a, a, a, b, b, a)
    weights = [sp.S.Zero]*6 + [sp.S.One]
    for i in range(5, -1, -1):
        weights[i] = reduce(-sum(weights[j] * columns[i][j-i] for j in range(i+1, 7)))
    assert all(reduce(sum(weights[j]*columns[i][j-i] for j in range(i,7))) == 0
               for i in range(6))
    alpha, beta = [reduce(sum(weights[j]*b[j-i] for j in range(i,7))) for i in (0,1)]
    assert beta != 0
    assert reduce(alpha + t*beta) == 0
    assert reduce(sum(weights)) == 0

    # Compare the two explicit completions of a zero baseline. In the all-one
    # case, the identity forces A1=t*A0 and the remaining rows are triangular.
    aa, bb = [u], [1-u]
    for k in range(1, 7):
        s = reduce(1-sum(aa[i]*b[k-i]+bb[i]*a[k-i] for i in range(k)))
        left = t*u if k == 1 else s if k in (3,4) else 0
        aa.append(left)
        bb.append(reduce(s-left))
    assert all(reduce(aa[k]-u*a[k]) == reduce(bb[k]-(1-u)*b[k]) == 0 for k in range(7))

    prime_data = []
    for modulus, image in ((163, -6), (1279, 168)):
        assert sp.isprime(modulus)
        assert all(int(f.subs(t, image)) % modulus == 0 for f in (h, alpha, beta))
        residues = [int(weight.subs(t, image)) % modulus for weight in weights]
        classes = defaultdict(list)
        for code in product((0,1), repeat=7):
            key = sum(x * y for x, y in zip(residues, code, strict=True)) % modulus
            classes[key].append("".join(map(str, code)))
        collisions = {key: words for key, words in classes.items() if len(words) > 1}
        prime_data.append({"prime": modulus, "t_image": image, "weights": residues,
                           "classes": len(classes), "nontrivial_collisions": collisions})
    assert prime_data[1]["nontrivial_collisions"] == {0: ["0000000", "1111111"]}
    result = {"status": "verified", "minimum_degree_for_identity": 15,
              "basis": "Z[t]/(H)", "lambda": [str(f) for f in weights], "alpha": str(alpha),
              "beta": str(beta), "binary_words_per_prime": 128, "prime_data": prime_data,
              "sympy": sp.__version__}
    with args.out.open("x") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
    print(json.dumps({"status": result["status"], "classes": [r["classes"] for r in prime_data]}))


if __name__ == "__main__":
    main()
