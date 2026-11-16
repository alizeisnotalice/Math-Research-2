"""Exact coordinate audit for H16(x^3)H16(x^4)."""
S = {0, 1, 2, 15, 16}


def main():
    target = {3*i + 4*j: (i, j) for i in S for j in S}
    assert len(target) == 25
    summed = {i + j for i in S for j in S}
    extras = []
    for n, base in target.items():
        for I in summed:
            for J in summed:
                if (I, J) != base and 3*I + 4*J == n:
                    extras.append((base, (I, J), (I-base[0])//4))
    assert len(extras) == 10
    assert all(t == 4 and j in {15, 16} and J == j-12
               for ((i, j), (I, J), t) in extras)
    print({"status": "passed_exact_H16_3_4_carry_audit", "carry_count": len(extras)})


if __name__ == "__main__":
    main()
