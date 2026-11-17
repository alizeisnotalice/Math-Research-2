#!/usr/bin/env python3
"""Strict finite-record audit of the completed remote r025 fixed-root jobs."""
import csv
import hashlib
import json
import re
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'remote-results'
SOURCE = ROOT / 'fixed_root_n16.cpp'
EXPECTED_SOURCE = 'd8e7040f251c60ba7b660fcdba9e153b9dc21e351105fb474ad32314cc3456a4'
TASKS = [(119, 262144, 128, 16, 12, 0xA024E1000000), (0, 262144, 32, 16, 3, 0xA024E2000000)]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def resource(path):
    text = path.read_text()
    fields = {}
    for label, key in [('real', 'wall_seconds'), ('maximum resident set size', 'maximum_RSS_bytes'), ('swaps', 'swaps'), ('peak memory footprint', 'peak_memory_footprint_bytes')]:
        match = re.search(r'^\s*([0-9.]+)\s+' + re.escape(label), text, re.MULTILINE)
        assert match, (path, label)
        fields[key] = match.group(1)
    return fields

def run(root, particles, reps, moves, threads, base):
    stem = f'n16-root{root}-p{particles}-r{reps}-l{moves}'
    csv_path = DATA / (stem + '.csv')
    metadata_path = DATA / (stem + '.csv.metadata.json')
    exit_path = DATA / (stem + '.exit.json')
    resource_path = DATA / (stem + '.resource.txt')
    metadata = json.loads(metadata_path.read_text())
    assert json.loads(exit_path.read_text())['exit_code'] == 0
    assert tuple(int(metadata[k]) for k in ('n', 'N', 'root', 'P', 'R', 'L', 'threads')) == (16, 256, root, particles, reps, moves, threads)
    assert int(metadata['seed_base']) == base
    with csv_path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == reps and {int(row['rep']) for row in rows} == set(range(reps))
    values, hashes = [], []
    with localcontext() as ctx:
        ctx.prec = 90
        for row in rows:
            rep = int(row['rep'])
            assert tuple(int(row[k]) for k in ('n', 'N', 'root', 'L')) == (16, 256, root, moves)
            assert int(row['seed']) == (base + root * 1000000007 + rep * 1009 + moves * 1000003) % (1 << 64)
            assert int(row['dead']) == 0
            assert int(row['proposals']) == particles * moves * 254
            assert 0 <= int(row['accept']) <= int(row['proposals'])
            h = int(row['terminal_hash'])
            assert 0 <= h < (1 << 64)
            hashes.append(h)
            z = Decimal(row['probZ'])
            assert z.is_finite() and z > 0
            values.append(z)
        total = sum(values)
        mean = total / reps
        delta = (sum((z - mean) ** 2 for z in values) / Decimal(reps * (reps - 1))).sqrt()
        return {
            'root': root, 'root_xy': [root % 16, root // 16],
            'metadata': metadata, 'exit': json.loads(exit_path.read_text()),
            'file_sha256': {p.name: sha(p) for p in (csv_path, metadata_path, exit_path, resource_path)},
            'resource': resource(resource_path),
            'mean_probZ_h_over_255_factorial': str(mean), 'replicate_delta_SE': str(delta),
            'replicate_weight_ESS': str(total * total / sum(v * v for v in values)),
            'max_replicate_weight_share': str(max(values) / total),
            'minimum_particle_parent_ESS': min(float(row['min_parent_ess']) for row in rows),
            'terminal_hashes': {'records_checked': len(hashes), 'all_uint64': True, 'distinct_count': len(set(hashes))},
        }

assert sha(SOURCE) == EXPECTED_SOURCE
driver = json.loads((DATA / 'driver.exit.json').read_text())
assert driver['exit_code'] == 0
out = {
    'status': 'passed',
    'scope': 'Finite SMC raw-record audit. The means and delta SEs are diagnostics, not confidence intervals or analytic inequalities.',
    'source_sha256': EXPECTED_SOURCE,
    'driver_exit': driver,
    'checks': ['driver and task exits zero', 'complete replicate sets', 'metadata and deterministic seed formulas', 'finite positive normalizers', 'dead=0', 'proposal/acceptance fields', 'terminal-hash uint64 field for every record'],
    'runs': [run(*task) for task in TASKS],
}
ROOT.joinpath('audit-remote-results.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
