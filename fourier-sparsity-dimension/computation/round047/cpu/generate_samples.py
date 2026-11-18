"""Generate 64 deterministic, independently sampled degree-three ANFs."""
import random
from pathlib import Path
import sys

seed = 202611180958
rng = random.Random(seed)
monomials = [a for a in range(1 << 9) if a.bit_count() <= 3]
assert len(monomials) == 130
rows = []
while len(rows) < 64:
    coefficients = [rng.getrandbits(1) for a in monomials]
    if not any(b and a.bit_count() == 3 for a, b in zip(monomials, coefficients)):
        continue
    rows.append("".join(map(str, coefficients)))
assert len(set(rows)) == 64
Path(sys.argv[1]).write_text("\n".join(rows) + "\n", encoding="ascii")
