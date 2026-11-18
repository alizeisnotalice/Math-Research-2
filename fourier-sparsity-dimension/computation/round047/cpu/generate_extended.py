"""Create disjoint 64-sample cohorts for the longer independent census."""
import json
from pathlib import Path
import random
import sys


directory = Path(sys.argv[1])
batches = int(sys.argv[2])
assert batches > 0
directory.mkdir()
seed = 202611181006
random_source = random.Random(seed)
terms = [x for x in range(512) if x.bit_count() <= 3]
seen = set()
for batch in range(batches):
    rows = []
    while len(rows) < 64:
        digits = [random_source.getrandbits(1) for _ in terms]
        if not any(value and term.bit_count() == 3 for term, value in zip(terms, digits)):
            continue
        row = "".join(map(str, digits))
        if row not in seen:
            seen.add(row)
            rows.append(row)
    (directory / f"cohort{batch:02d}.txt").write_text("\n".join(rows) + "\n", encoding="ascii")
(directory / "manifest.json").write_text(
    json.dumps({"seed": seed, "batches": batches, "samples": len(seen)}, indent=2) + "\n",
    encoding="ascii",
)
