#!/usr/bin/env python3
"""Portable finite-ledger audit for the three archived n=11 root runs."""
import csv
import hashlib
import json
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
N, P, R, L = 121, 262144, 128, 16
SEED0 = 0xA024B100
ROOTS = {13: (2, 1), 24: (2, 2), 25: (3, 2)}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def one(root, xy):
    stem = HERE / f'n11-root{root}-p{P}-r{R}-l{L}'
    csv_path = stem.with_suffix('.csv')
    meta = json.loads(Path(str(csv_path) + '.metadata.json').read_text())
    ended = json.loads(stem.with_suffix('.exit.json').read_text())
    assert ended['exit_code'] == 0
    assert tuple(int(meta[k]) for k in ('n', 'N', 'root', 'P', 'R', 'L')) == (11, N, root, P, R, L)
    assert tuple(meta['root_xy']) == xy and int(meta['threads']) == 5
    base = SEED0 + root * (1 << 36)
    assert int(meta['seed_base']) == base
    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == R and {int(row['rep']) for row in rows} == set(range(R))
    with localcontext() as ctx:
        ctx.prec = 90
        vals = []
        for row in rows:
            rep = int(row['rep'])
            assert int(row['seed']) == (base + root * 1000000007 + rep * 1009 + L * 1000003) % (1 << 64)
            assert int(row['dead']) == 0
            assert int(row['proposals']) == P * L * (N - 2)
            assert 0 <= int(row['accept']) <= int(row['proposals'])
            z = Decimal(row['Z'])
            assert z.is_finite() and z > 0
            vals.append(z)
        total = sum(vals)
        mean = total / R
        se = (sum((z - mean) ** 2 for z in vals) / Decimal(R * (R - 1))).sqrt()
        return {'root': root, 'xy': xy, 'csv_sha256': sha(csv_path),
                'metadata_sha256': sha(Path(str(csv_path) + '.metadata.json')),
                'exit_sha256': sha(stem.with_suffix('.exit.json')),
                'resource_sha256': sha(stem.with_suffix('.resource.txt')),
                'mean_h_root': str(mean), 'replicate_delta_SE': str(se),
                'replicate_weight_ESS': str(total * total / sum(z*z for z in vals)),
                'max_replicate_weight_share': str(max(vals) / total),
                'minimum_particle_parent_ESS': min(float(row['min_parent_ess']) for row in rows),
                'accepted': sum(int(row['accept']) for row in rows),
                'proposals': sum(int(row['proposals']) for row in rows)}

out = {'scope': 'Finite record audit only; no root-transport or confidence claim.',
       'checks': ['exit code zero', '128 distinct replicates', 'metadata',
                  'deterministic seeds', 'positive finite Z', 'dead=0',
                  'exact proposal count', 'acceptance range'],
       'runs': [one(root, xy) for root, xy in ROOTS.items()]}
(HERE / 'verify-archive.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
