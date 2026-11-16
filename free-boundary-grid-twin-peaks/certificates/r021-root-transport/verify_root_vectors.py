#!/usr/bin/env python3
"""Check exact root vectors emitted by the independent-set DP for n=2,...,6."""
from fractions import Fraction
import json
from pathlib import Path


BASE = Path(__file__).parent
EXPECTED = {
    2: Fraction(1),
    3: Fraction(11, 4),
    4: Fraction(19681, 3560),
    5: Fraction(6786834357579, 719889959846),
    6: Fraction(68708495735539617742323275, 4793987217870606674220112),
}


def main():
    out = []
    for n, expected in EXPECTED.items():
        data = json.loads((BASE / f"n{n}-roots.json").read_text())
        roots = list(map(int, data["single_peak_root_counts"]))
        assert len(roots) == n * n
        assert sum(roots) == int(data["A1"])
        # D4 symmetry is an independent structural check on the root vector.
        at = lambda x, y: roots[y * n + x]
        for y in range(n):
            for x in range(n):
                assert at(x, y) == at(n - 1 - x, y) == at(x, n - 1 - y)
        inward_edges = 0
        for y in range(n):
            for x in range((n - 1) // 2):
                assert at(x + 1, y) >= at(x, y)
                inward_edges += 1
        for x in range(n):
            for y in range((n - 1) // 2):
                assert at(x, y + 1) >= at(x, y)
                inward_edges += 1
        max_ratio = max(
            Fraction(max(at(x, y), at(xx, yy)), min(at(x, y), at(xx, yy)))
            for y in range(n) for x in range(n)
            for xx, yy in ((x + 1, y), (x, y + 1))
            if xx < n and yy < n
        )
        assert max_ratio == expected
        # This is only a finite audit of the proposed all-row transport scale
        # a[x+1,y]/a[x,y] <= C*n^2/(x+1)^2.  Record the least C on this
        # finite grid rather than treating it as an asymptotic statement.
        scaled = [
            Fraction(at(x + 1, y), at(x, y)) * Fraction((x + 1) ** 2, n * n)
            for y in range(n) for x in range(n - 1)
        ]
        transport_constant = max(scaled)
        assert transport_constant <= 1
        binomial_scaled = [
            Fraction(at(x + 1, y), at(x, y))
            / Fraction((n - 1 - x) ** 2, (x + 1) ** 2)
            for y in range(n) for x in range((n - 1) // 2)
        ]
        # Again this is an exact finite check only.  The candidate is
        # a[x+1,y]/a[x,y] <= ((n-1-x)/(x+1))^2 before the central line.
        binomial_scale_constant = max(binomial_scaled, default=Fraction(0))
        assert binomial_scale_constant <= 1
        central_edges = 0
        if n % 2 == 0:
            x = n // 2 - 1
            for y in range(n):
                assert at(x + 1, y) == at(x, y)
                central_edges += 1
        out.append({"n": n, "A1": str(sum(roots)), "max_adjacent_ratio": str(max_ratio),
                    "inward_coordinate_edges_checked": inward_edges,
                    "all_row_transport_scale_constant": str(transport_constant),
                    "strict_inward_binomial_square_scale_constant": str(binomial_scale_constant),
                    "central_horizontal_reflection_edges_checked": central_edges})
    payload = {"status": "passed", "runs": out}
    (BASE / "root-vector-verified.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
