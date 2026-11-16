#!/usr/bin/env python3
"""Recheck the retained completed n=16 fixed-root production CSV."""
import csv
import hashlib
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
RAW = HERE / 'raw' / 'n16-root119-p65536-r64-l16-direct.csv'
META = pathlib.Path(str(RAW) + '.metadata.json')
HASHES = HERE / 'source' / 'production-direct.sha256'

expected = {}
for line in HASHES.read_text().splitlines():
    digest, name = line.split(maxsplit=1)
    expected[name] = digest
assert hashlib.sha256((HERE / 'source' / 'fixed_root_n16.cpp').read_bytes()).hexdigest() == expected['fixed_root_n16.cpp']
assert hashlib.sha256(RAW.read_bytes()).hexdigest() == expected[RAW.name]
assert hashlib.sha256(META.read_bytes()).hexdigest() == expected[META.name]
m = json.loads(META.read_text())
assert (m['n'], m['N'], m['root'], m['P'], m['R'], m['L'], m['threads']) == (16, 256, 119, 65536, 64, 16, 12)
rows = list(csv.DictReader(RAW.open()))
assert len(rows) == 64 and {int(x['rep']) for x in rows} == set(range(64))
base = int(m['seed_base'])
for x in rows:
    r = int(x['rep'])
    assert int(x['seed']) == (base + 119 * 1000000007 + r * 1009 + 16 * 1000003) % (1 << 64)
    assert int(x['dead']) == 0 and float(x['probZ']) > 0 and float(x['min_parent_ess']) > 0
    assert int(x['proposals']) == 65536 * 16 * 254
print(json.dumps({'passed': True, 'raw_rows': len(rows),
                  'source_sha256': expected['fixed_root_n16.cpp'],
                  'raw_sha256': expected[RAW.name]}, indent=2))
