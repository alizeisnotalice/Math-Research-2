"""Recover every rational coordinate of the eight witnesses and the optional ninth set."""
import argparse
import gzip
import hashlib
import json
from itertools import accumulate
from pathlib import Path

EXPECTED = "f1d72553f8fcb40f89ac7c85cd09c2a4a81ecc0ddf380d1da23961f17cbe4f10"


def restore(values):
    for _ in range(4):
        values = list(accumulate(values))
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    source = Path(__file__).with_name("points-d4.json.gz")
    with gzip.open(source, "rt", encoding="utf-8") as stream:
        records = json.load(stream)
    point_sets = []
    for body, denominator, x_differences, y_differences in records:
        vertices = list(zip(restore(x_differences), restore(y_differences), strict=True))
        point_sets.append(dict(body=body, denominator=denominator, points=vertices))
    raw = (json.dumps(point_sets, separators=(",", ":")) + "\n").encode()
    if hashlib.sha256(raw).hexdigest() != EXPECTED:
        raise ValueError("Decoded coordinate table differs from the verified witness.")
    args.output.write_bytes(raw)
    print(EXPECTED, [len(group["points"]) for group in point_sets])


if __name__ == "__main__":
    main()
