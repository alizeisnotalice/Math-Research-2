"""Enumerate every intrinsic parity direction of the first aperiodic certificate."""
import json
import sys
from collections import Counter
from pathlib import Path

from aperiodic_graph_audit import rank


def main(reference):
    row = next(row for row in json.loads(reference.read_text())["rows"] if row["periods"] == [0])
    origin = row["support"][0]
    spectrum = {k ^ origin: c for k, c in zip(row["support"], row["coefficient_numerators_over_8"])}
    space = {0}
    for v in spectrum:
        if v not in space:
            space.update([u ^ v for u in space])
    assert len(space) == 2**14
    outcomes = Counter()
    for query in space - {0}:
        projected = {min(v, v ^ query) for v in spectrum}
        if len(projected) == len(spectrum):
            # Both equal supports have the same augmented dimension.
            h = rank(v | 32768 for v in projected)
            assert h == 14
            s0 = s1 = len(projected)
        else:
            left, right = [], []
            for v in projected:
                a, b = spectrum.get(v, 0), spectrum.get(v ^ query, 0)
                if a + b:
                    left.append(v | 32768)
                if a - b:
                    right.append(v | 32768)
            h = rank(left) + rank(right) - rank(left + right)
            s0, s1 = len(left), len(right)
        assert s0 == s1 and (16 + h)**2 < 16*s0
        outcomes[f"{s0},{h}"] += 1
    print(json.dumps(dict(sorted(outcomes.items())), indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
