"""Exact integer carry ledger for H10(x^3)H10(x^4)."""
from collections import defaultdict
import json

S = (0, 1, 2, 9, 10)
A, B = 3, 4


def main():
    canonical = {(A * i + B * j): (i, j) for i in S for j in S}
    assert len(canonical) == 25
    summands = {p + q for p in S for q in S}

    extras = []
    for n, base in sorted(canonical.items()):
        for I in sorted(summands):
            for J in sorted(summands):
                if (I, J) != base and A * I + B * J == n:
                    t = (I - base[0]) // B
                    assert I == base[0] + B * t and J == base[1] - A * t
                    extras.append({"canonical": list(base), "noncanonical": [I, J], "t": t})

    expected = [
        {"canonical": [0, 9], "noncanonical": [12, 0], "t": 3},
        {"canonical": [1, 9], "noncanonical": [9, 3], "t": 2},
        {"canonical": [0, 10], "noncanonical": [12, 1], "t": 3},
        {"canonical": [2, 9], "noncanonical": [10, 3], "t": 2},
        {"canonical": [1, 10], "noncanonical": [9, 4], "t": 2},
        {"canonical": [2, 10], "noncanonical": [10, 4], "t": 2},
        {"canonical": [10, 9], "noncanonical": [18, 3], "t": 2},
        {"canonical": [10, 10], "noncanonical": [18, 4], "t": 2},
    ]
    assert extras == expected, extras

    print(json.dumps({
        "status": "passed_exact_H10_3_4_multicarry_audit",
        "S": list(S),
        "target_terms": len(canonical),
        "target_binary": True,
        "actual_product_coordinate_set": sorted(summands),
        "noncanonical_carries": extras,
        "carry_t_values": sorted({row["t"] for row in extras}),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
