#!/usr/bin/env python3
"""Exact four-variable check of extensions of a fixed completed E4."""

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
import time


def walsh(mask, support):
    points = [x for x in range(16) if (support >> x) & 1]
    return tuple(
        sum(1 - 2 * ((((mask >> x) & 1) ^ ((x & frequency).bit_count() & 1)))
            for x in points)
        for frequency in range(16)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()

    characters = tuple(
        sum(((x & frequency).bit_count() & 1) << x for x in range(16))
        for frequency in range(16)
    )
    bent = {
        phase for phase in range(65536)
        if all((phase ^ character).bit_count() in (6, 10) for character in characters)
    }
    base = sum(((((x >> 0) & 1) * ((x >> 1) & 1))
                ^ (((x >> 2) & 1) * ((x >> 3) & 1))) << x for x in range(16))
    flat = 255  # x3=0, common half-block direction x2.
    e4 = {(characters[direction] & flat) ^ offset
          for direction in (0, 1, 2, 3) for offset in (0, flat)}
    assert all(base ^ word in bent for word in e4)
    compatible = [direction for direction in range(65536)
                  if all(base ^ word ^ direction in bent for word in e4)]
    old_support = {frequency for frequency, value in enumerate(walsh(base, flat)) if value}
    branch_counts = Counter()
    weight_counts = {0: Counter(), 1: Counter()}

    for direction in compatible:
        if direction in e4:
            continue
        on_flat = direction & flat
        outside = direction ^ on_flat
        slope = (on_flat ^ (on_flat >> 4)) & 1
        assert on_flat in {character & flat for character in characters} | {
            (character & flat) ^ flat for character in characters
        }
        assert outside.bit_count() >= 4
        branch_counts[slope] += 1
        weight_counts[slope][outside.bit_count()] += 1
        new_support = {frequency for frequency, value
                       in enumerate(walsh(base ^ direction, flat)) if value}
        if slope == 0:
            assert old_support == new_support and base ^ outside in bent
        else:
            assert not old_support & new_support
            changed = old_support | new_support
            outside_spectrum = walsh(base, outside)
            assert len(changed) == 16 and all(abs(outside_spectrum[freq]) == 2
                                              for freq in changed)
            assert outside.bit_count() == 4
            points = [point for point in range(16) if (outside >> point) & 1]
            assert points[0] ^ points[1] ^ points[2] ^ points[3] == 0

    disjoint_q = sum(1 << point for point in range(16) if point & 9 == 8)
    assert disjoint_q.bit_count() == 4 and disjoint_q in compatible
    assert min(weight_counts[0]) == min(weight_counts[1]) == 4
    result = {
        "ambient_dimension": 4,
        "q": 4,
        "base_phase": "x0*x1 + x2*x3",
        "flat": "x3=0",
        "total_directions": 65536,
        "bent_phases": len(bent),
        "compatible_directions_including_e4": len(compatible),
        "compatible_extensions_by_slope": dict(sorted(branch_counts.items())),
        "outside_weight_by_slope": {slope: dict(sorted(counts.items()))
                                     for slope, counts in weight_counts.items()},
        "disjoint_q_mask": disjoint_q,
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers, no randomness",
        "elapsed_seconds": time.perf_counter() - started,
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"checked all {result['total_directions']} directions; {len(compatible)} compatible")


if __name__ == "__main__":
    main()
