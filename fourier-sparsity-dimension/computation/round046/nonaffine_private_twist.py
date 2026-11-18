#!/usr/bin/env python3
"""Verify a counterexample with one nonlinear private shift."""

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import runpy
import sys


source = Path(sys.argv[1])
fwht = runpy.run_path(str(source))["fwht"]


def route(x, w, width):
    nonlinear_private_root = (x & 1) ^ (w != 0)
    shared_is_in_E = all((w >> j) & 1 == 0 for j in range(width))
    if (nonlinear_private_root == 0) == shared_is_in_E:
        return 1 - 2 * ((x >> 1) & 1), nonlinear_private_root
    return 1, nonlinear_private_root


def example(width):
    W = 1 << width
    truth = [route(x, w, width)[0] for x in range(4) for w in range(W)]
    assert truth == [1 if x & 1 else 1 - 2 * (x >> 1) for x in range(4) for _ in range(W)]
    support = [i for i, coefficient in enumerate(fwht(truth)) if coefficient]
    assert support == [0, W, 2 * W, 3 * W]
    band_support = [i for i in support if i // W in (2, 3)]
    assert len(band_support) == 2

    root_deltas = {}
    for w in range(W):
        conditional = [
            [out for x in range(4) for out, branch in [route(x, w, width)] if branch == bit]
            for bit in (0, 1)
        ]
        assert [len(c) for c in conditional] == [2, 2]
        root_deltas[w] = Fraction(sum(conditional[0]) - sum(conditional[1]), 4)
    assert set(root_deltas.values()) == {Fraction(-1, 2), Fraction(1, 2)}
    assert route(0, 0, width) == (1, 0) and route(2, 0, width) == (-1, 0)

    # Left-side departures from E have constant leaves. The right side contains
    # sensitive x_1 queries, so the root is the only blocker for the left target.
    for w in range(1, W):
        assert route(1, w, width) == route(3, w, width) == (1, 0)
        assert route(0, w, width) == (1, 1) and route(2, w, width) == (-1, 1)
    p = Fraction(1, W)
    local_false_lower_bound = 1 / p  # k=ell=1.
    two_terms = 1 + Fraction(1, 2) / p  # Root has k=ell=0 and is marked.
    assert local_false_lower_bound > len(band_support)
    return {
        "shared_bits": width, "full_Walsh_support": support, "sparsity": len(support),
        "target_private_depth": 1, "target_ell": 1, "target_p": str(p),
        "target_band_support": band_support, "claimed_local_lower_bound": str(local_false_lower_bound),
        "root_conditional_half_differences": {w: str(d) for w, d in root_deltas.items()},
        "claimed_global_bound_from_only_two_terms": str(two_terms),
        "already_violates_global_bound": two_terms > len(support),
    }


result = {
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "round045_helper_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "cases": [example(width) for width in (2, 3)],
}
assert result["cases"][1]["already_violates_global_bound"]
Path(sys.argv[2]).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
