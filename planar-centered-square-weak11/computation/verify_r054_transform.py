#!/usr/bin/env python3
"""Verify a saved left-reflected input pointwise with exact rational arithmetic."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

p = argparse.ArgumentParser()
p.add_argument('source', type=Path)
p.add_argument('reflected', type=Path)
a = p.parse_args()
s, t = (json.loads(f.read_text()) for f in (a.source, a.reflected))
assert set(s) == set(t) == {'lam', 'Px', 'Py', 'pts'}
assert all(F(s[k]) == F(t[k]) for k in ('lam', 'Px', 'Py'))
assert len(s['pts']) == len(t['pts']) == 204
X = max(F(pt[0]) for pt in s['pts'])
assert all((X-F(x), F(y), F(m)) == tuple(map(F, q)) for (x, y, m), q in zip(s['pts'], t['pts']))
mass = sum((F(pt[2]) for pt in s['pts']), F(0))
assert mass == sum((F(pt[2]) for pt in t['pts']), F(0))
print(json.dumps({'status':'PASS', 'python':sys.version, 'source_sha256':hashlib.sha256(a.source.read_bytes()).hexdigest(), 'reflected_sha256':hashlib.sha256(a.reflected.read_bytes()).hexdigest(), 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'points':len(s['pts']), 'X':str(X), 'mass':str(mass), 'scope':'pointwise exact transform only; no seam, halo, or weak-type bound verified'}, indent=2))
