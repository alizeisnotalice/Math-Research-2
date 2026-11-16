#!/usr/bin/env python3
"""Exact arithmetic audit for p2q7_explicit_certificate.md."""


def main() -> None:
    p, q, r = 2, 7, 10_000
    h_indices = {r - 2, r - 6, r - 7, r - 21}
    assert all(x < r for x in h_indices)
    assert (r - 2) - (r - 6) == p * 2
    assert (r - 7) - (r - 21) == q * 2
    assert p * (r - 7) + q * (r - 6) == p * (r - 21) + q * (r - 2)
    print("four-offset C-collision certificate passed")


if __name__ == "__main__":
    main()
