"""Exact finite audit of the N^2 uniqueness lemma on representative inputs."""


def all_pairs(n: int, a: int, b: int):
    return [(i, j) for i in range(n // a + 1) for j in range(n // b + 1)
            if a * i + b * j == n]


def check(a: int, b: int, ip, jq):
    target = {a * i + b * j for i in ip for j in jq}
    assert len(target) == len(ip) * len(jq)
    for n in target:
        pairs = all_pairs(n, a, b)
        assert len(pairs) == 1, (a, b, n, pairs)


def main():
    # H support with the smallest permitted coprime pair and a second pair.
    H = {0, 1, 2, 5, 6}
    check(7, 8, H, H)
    check(8, 9, H, H)
    # Boundary a=1 forces the second mother to be constant.
    check(1, 3, {0, 1, 2}, {0})
    print({"status": "passed", "examples": ["H_a7_b8", "H_a8_b9", "a1_Qconstant"]})


if __name__ == "__main__":
    main()
