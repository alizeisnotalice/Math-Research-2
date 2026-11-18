#!/usr/bin/env python3
"""Check every four-variable phase and all cuts on one size-eight flat."""

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
import time


def character_masks(size):
    return tuple(
        sum(((point & frequency).bit_count() & 1) << point for point in range(size))
        for frequency in range(size)
    )


def support(mask, characters, size, amplitude):
    spectrum = tuple(size - 2 * (mask ^ character).bit_count() for character in characters)
    assert all(value in (-amplitude, 0, amplitude) for value in spectrum)
    return {index for index, value in enumerate(spectrum) if value}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    full = character_masks(16)
    flat = character_masks(8)
    bent = {
        mask
        for mask in range(1 << 16)
        if all((mask ^ character).bit_count() in (6, 10) for character in full)
    }
    flat_indicator = 255
    bases = sorted(mask for mask in bent if mask ^ flat_indicator in bent)
    direction_histogram = Counter()
    legal_count = 0
    for base in bases:
        old_support = support(base & flat_indicator, flat, 8, 4)
        origin = min(old_support)
        directions = {frequency ^ origin for frequency in old_support}
        assert len(directions) == 4
        assert all(left ^ right in directions for left in directions for right in directions)
        predicted = {flat[form] ^ constant for form in directions for constant in (0, flat_indicator)}
        actual = {cut for cut in range(256) if base ^ cut in bent}
        assert len(actual) == 8 and actual == predicted
        for cut in actual:
            assert support((base ^ cut) & flat_indicator, flat, 8, 4) == old_support
        assert Counter(cut.bit_count() for cut in actual) == {0: 1, 4: 6, 8: 1}
        legal_count += len(actual)
        direction_histogram[",".join(map(str, sorted(directions - {0})))] += 1
    result = {
        "ambient_dimension": 4,
        "flat": "x3=0; x0,x1,x2 free",
        "total_phases": 1 << 16,
        "bent_phases": len(bent),
        "bases_with_legal_full_flip": len(bases),
        "tested_cut_pairs": len(bases) * 256,
        "legal_cut_pairs": legal_count,
        "direction_histogram": dict(sorted(direction_histogram.items())),
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers, no randomness",
        "elapsed_seconds": time.perf_counter() - started,
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    assert len(bent) == len(bases) == 896
    assert len(direction_histogram) == 7 and set(direction_histogram.values()) == {128}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"verified {result['tested_cut_pairs']} base-cut pairs")


if __name__ == "__main__":
    main()
