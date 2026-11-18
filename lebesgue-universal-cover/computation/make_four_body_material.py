"""Generate exact, strictly interior integer material at each nominal four-body angle."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json
import numpy as np
from flint import arb, ctx

COORD = 2**28
FINE = 2**44
GRID = 1024


def rational(x):
    z = Q(x)
    return arb(z.numerator) / z.denominator


def make(header, subdiv, disk):
    ctx.prec = 128
    hi = header["hi"]
    counts = [disk, 3 * subdiv, 5 * subdiv, 7 * subdiv, 0]
    table = np.empty((GRID + 1, sum(counts), 2), dtype=np.int32)
    offset = 0
    for body, n in enumerate((0, 3, 5, 7)):
        total = counts[body]
        radius = 1 / (2 * (arb.pi() / (2 * n)).cos()) if n else arb(1) / 2
        corners = []
        for j in range(n):
            s, c = (2 * arb.pi() * j / n).sin_cos()
            corners.append((radius * c, radius * s))
        base = []
        for j in range(total):
            s, c = (2 * arb.pi() * j / total).sin_cos()
            if n:
                m = (Q(j * n, total) + Q(n, 2) + Q(1, 2)).__floor__() % n
                v = c * corners[m][0] + s * corners[m][1]
                rad = v + (1 - radius**2 + v**2).sqrt()
            else:
                rad = radius
            rad *= 1 - arb(1) / 10**6
            base.append((rad * c, rad * s))
        cap = rational(hi[3 * body - 4]) if body >= 2 else arb(0)
        for i in range(GRID + 1) if body >= 2 else [0]:
            s, c = (cap * i / GRID).sin_cos()

            def rotate(p):
                x, y = p
                return c * x - s * y, s * x + c * y

            vv = list(map(rotate, corners))
            for j, p in enumerate(base):
                x, y = rotate(p)
                a, b = (round(float(z.mid()) * COORD) for z in (x, y))
                X, Y = arb(a) / COORD, arb(b) / COORD
                if not n:
                    assert X * X + Y * Y < arb(1) / 4
                else:
                    assert all((X - u) ** 2 + (Y - v) ** 2 < 1 for u, v in vv)
                table[i, offset + j] = a, b
        if body < 2:
            table[1:, offset : offset + total] = table[0, offset : offset + total]
        offset += total
    assert np.max(np.abs(table.astype(np.int64))) < COORD
    # Last three entries only maintain a fixed-width native interface.
    original = list(hi) + [1 / 2, 0, 0]
    assert all(Q(v) * FINE == int(Q(v) * FINE) for v in original)
    return dict(
        points=table,
        counts=np.asarray(counts, dtype=np.int32),
        original_hi=np.array([int(Q(v) * FINE) for v in original], dtype=np.int64),
    )


if __name__ == "__main__":
    ap = ArgumentParser()
    ap.add_argument("tree", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--subdiv", type=int, default=64)
    ap.add_argument("--disk", type=int, default=384)
    a = ap.parse_args()
    h = json.loads((a.tree / "part0000.json").read_bytes())
    assert h["ns"] == [3, 5, 7] and len(h["hi"]) == 8
    data = make(h, a.subdiv, a.disk)
    np.savez_compressed(a.output, **data)
    print(
        json.dumps(
            dict(
                counts=data["counts"].tolist(),
                points_sha256=sha256(data["points"].tobytes()).hexdigest(),
                file_sha256=sha256(a.output.read_bytes()).hexdigest(),
            )
        )
    )
