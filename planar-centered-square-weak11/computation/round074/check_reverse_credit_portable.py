#!/usr/bin/env python3
"""Portable exact equality check from the public compact R74 summary only.

Adapted from the separately archived independent diagonal checker. This
replays one seed and exact three-subset geometry, not the 1200-case search.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import platform
import random


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--summary', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    summary = json.loads(args.summary.read_text())
    top = summary['top20'][0]
    assert (top['count'], top['seed'], top['mode']) == (2, 742090, 1)
    rng = random.Random(top['seed'])
    pts = set()
    while len(pts) < 2:
        x, y = rng.choice((-10, 10)), rng.randint(-10, 10)
        if rng.randrange(2):
            x, y = y, x
        if (x, y) != (0, 0):
            pts.add((x, y))
    pts = sorted(pts)
    weights = [rng.randint(1, 20) for _ in pts]
    assert weights == [8, 7]
    points = [tuple(F(t, 20) for t in p) for p in pts]
    assert points == [(F(-1, 2), F(-1, 2)), (F(1, 2), F(1, 2))]
    scale = F(1, sum(weights))
    masses = [scale * w for w in weights]
    assert [list(map(str, p)) for p in points] == top['positions']
    assert list(map(str, masses)) == top['masses']
    # All points have distance 1/2, so this is the entire tied prefix menu.
    assert all(max(map(abs, p)) == F(1, 2) for p in points)
    assert sum(masses) == 1 == 4 * F(1, 2)**2
    subsets = []
    for mask in (1, 2, 3):
        ids = [i for i in range(2) if (mask >> i) & 1]
        mass = sum(masses[i] for i in ids)
        span = max(max(points[i][j] for i in ids)
                   - min(points[i][j] for i in ids) for j in (0, 1))
        subsets.append({'mask': mask, 'mass': str(mass),
                        'span_squared': str(span**2),
                        'eligible': mass > span**2})
    assert [s['eligible'] for s in subsets] == [True, True, False]
    # Each halo radius is <1/2, hence the opposite corner squares are disjoint.
    # Intersecting each with closed Q leaves its quarter, of area mass/4.
    quarter_areas = [m / 4 for m in masses]
    assert all(0 < a < F(1, 4) for a in quarter_areas)
    area = sum(quarter_areas)
    assert area == F(1, 4) and sum(masses) / area == 4
    out = {'schema': 'r074-portable-exact-equality-v1', 'passed': True,
           'python': platform.python_version(), 'pid': os.getpid(),
           'summary_sha256': hashlib.sha256(args.summary.read_bytes()).hexdigest(),
           'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'seed': top['seed'], 'raw_weights': weights, 'scale': str(scale),
           'subsets': subsets, 'closed_Q_mass': '1', 'clipped_area': str(area),
           'exact_ratio': '4',
           'scope': 'One seed and all three subsets. No raw-output dependency; no other sample area replay and no general theorem.'}
    args.output.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'passed': True, 'exact_ratio': '4', 'pid': os.getpid()}))


if __name__ == '__main__':
    main()
