#!/usr/bin/env python3
"""Integer Walsh checks for two genuinely nonparallel affine partitions."""
import hashlib
import json
import platform
import time
from pathlib import Path
from verify_rank_two_families import basis, span


def truth(fn):
    return sum(fn(tuple((x >> j) & 1 for j in range(6))) << x for x in range(64))


def fwht(values):
    a = list(values)
    step = 1
    while step < len(a):
        for i in range(0, len(a), step * 2):
            for j in range(i, i + step):
                p, q = a[j], a[j + step]
                a[j], a[j + step] = p + q, p - q
        step *= 2
    return a


def case(name, coupled, mixed):
    # Domain bits: x,y,z,t,w,u. Base phase xy+zw+tu.
    f0 = truth(lambda b: (b[0] * b[1]) ^ (b[2] * b[4]) ^ (b[3] * b[5]))
    blocks = [[x for x in range(64) if ((x >> 4) & 1) == 0 and ((x >> 5) & 1) == c]
              for c in (0, 1)]
    blocks += [[x for x in range(64) if ((x >> 4) & 1) == 1 and ((x >> 3) & 1) == c]
               for c in (0, 1)]
    if mixed:
        a = blocks.pop(0)
        blocks += [[x for x in a if (x & 1) == c] for c in (0, 1)]
    assert sorted(x for a in blocks for x in a) == list(range(64))
    descriptions = []
    for a in blocks:
        p = min(a)
        h = span(x ^ p for x in a)
        assert h == frozenset(x ^ p for x in a)
        descriptions.append((p, basis(h), h))
    assert len({h for _, _, h in descriptions}) > 1
    indicators = [sum(1 << x for x in a) for a in blocks]
    directions = [truth(lambda b: 1)] + [truth(lambda b, j=j: b[j]) for j in range(6)] + indicators
    if coupled:
        directions += [truth(lambda b: b[0] * b[2]), truth(lambda b: b[1] * b[2])]
    bs = basis(directions)
    phases = sorted(f0 ^ v for v in span(bs))
    max_quad_rank = 0
    branch_checks = 0
    for phase in phases:
        w = fwht([1 - 2 * ((phase >> x) & 1) for x in range(64)])
        assert all(abs(c) == 8 for c in w)
        occupied = set()
        for a, (p, hb, h) in zip(blocks, descriptions):
            padded = fwht([1 - 2 * ((phase >> x) & 1) if x in h_p else 0
                           for h_p in [set(a)] for x in range(64)])
            supp = {i for i, val in enumerate(padded) if val}
            assert not occupied & supp
            occupied |= supp
            assert all(abs(padded[i]) == 8 for i in supp)
            points = [p]
            for d in hb:
                points += [x ^ d for x in points]
            intrinsic = fwht([1 - 2 * ((phase >> x) & 1) for x in points])
            assert sum(c != 0 for c in intrinsic) == (len(a) // 8) ** 2
            # Algebraic normal form and alternating-form rank on intrinsic coordinates.
            anf = [(phase >> x) & 1 for x in points]
            for j in range(len(hb)):
                for x in range(len(anf)):
                    if x & (1 << j):
                        anf[x] ^= anf[x ^ (1 << j)]
            assert all(not c or x.bit_count() <= 2 for x, c in enumerate(anf))
            rows = [0] * len(hb)
            for i in range(len(hb)):
                for j in range(i + 1, len(hb)):
                    if anf[(1 << i) | (1 << j)]:
                        rows[i] ^= 1 << j
                        rows[j] ^= 1 << i
            rank = len(basis(rows))
            assert rank == (2 if len(a) == 16 else 0)
            max_quad_rank = max(max_quad_rank, rank)
            branch_checks += 1
        assert occupied == set(range(64))
    return {"name": name, "bits": "x,y,z,t,w,u", "f0_mask": f0,
            "basis_masks": list(bs), "dimension": len(bs), "phase_count": len(phases),
            "blocks": blocks, "distinct_block_directions": len({h for _, _, h in descriptions}),
            "branch_checks": branch_checks, "max_intrinsic_quadratic_rank": max_quad_rank,
            "full_coefficients_checked": len(phases) * 64,
            "all_coefficients_exact_abs8": True}


def main():
    start = time.monotonic()
    cases = [case("nonparallel_coupled_four_double_blocks", True, False),
             case("nonparallel_mixed_three_double_two_minimal", False, True)]
    root = Path(__file__).parent
    result = {"schema": "a28-r017-nonparallel-exact-v1", "cases": cases,
              "seconds": time.monotonic() - start, "python": platform.python_version(),
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "helper_sha256": hashlib.sha256((root / "verify_rank_two_families.py").read_bytes()).hexdigest(),
              "arithmetic": "integer Walsh and exact GF(2)", "failed_assertions": 0}
    (root / "nonparallel-cases.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "cases"}))
    for c in cases:
        print(c["name"], c["dimension"], c["phase_count"], c["branch_checks"])


if __name__ == "__main__":
    main()
