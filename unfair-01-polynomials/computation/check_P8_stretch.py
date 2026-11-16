"""Exact integer audit for the P8 stretch proof; no numerical search."""
from __future__ import annotations

S = (0, 2, 3, 4, 6, 7, 8)


def support_c(L: int) -> set[int]:
    terms = [i + L * j for i in S for j in S]
    assert len(terms) == len(set(terms)) == 49
    return set(terms)


def main() -> None:
    for L in range(9, 201):
        c = support_c(L)
        assert 9 not in c and 10 not in c and 12 not in c
        for h in range(1, L):
            assert {(s + h) % L for s in S} != set(S)
    # Exact low-end convolution terms after the orientation a2=a3=1,
    # b2=b3=b6=b7=b8=0, a4=1-t, b4=t, a6=a7=1-t, a8=1.
    # c9=b6+b7, c10=b8+t(1-t), and at t=1 c12 includes a8*b4=1.
    print({"status": "passed", "L_checked": [9, 200],
           "support_size": 49, "zero_coefficients": [9, 10, 12],
           "nonzero_translation_stabilizers_found": 0})


if __name__ == "__main__":
    main()
