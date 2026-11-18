"""Verify the n=7 count ledger and rerun small and selected n=7 checks."""
import hashlib
import itertools
import json
import math
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

cert = Path(__file__).resolve().parent
project = cert.parents[1]
out = Path(sys.argv[1])
out.mkdir()
records = []

def run(command):
    process = subprocess.Popen(list(map(str, command)), stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    if process.returncode:
        raise RuntimeError((command, process.returncode, stdout, stderr))
    records.append(dict(pid=process.pid, returncode=0, waited_reaped=True))
    return stdout

norm = project / 'computation/results/n7-exact-r005.json'
a2 = int(json.loads(norm.read_text())['A2'])
assert hashlib.sha256(norm.read_bytes()).hexdigest() ==     'ff7974ff8487a339516fee95324674efd873abd9a36b8fc50b0f8121c632f8d6'
rows = {r['curve']: r for r in json.loads((cert / 'n7-curves.json').read_text())}
profiles = {''.join(chars) for chars in itertools.product('34', repeat=7)}
assert set(rows) == profiles
for h, row in rows.items():
    k = ''.join(str(7-int(c)) for c in h)
    assert row['complement_curve'] == k
    assert row['left_vertices'] == sum(map(int, h))
    assert row['right_vertices'] == sum(map(int, k))
    assert row['left_rooted_count'] == rows[k]['right_rooted_count']
    assert row['right_rooted_count'] == rows[k]['left_rooted_count']
    assert row['left_rooted_count'] == rows[h[::-1]]['left_rooted_count']
    expected = (math.comb(49, row['left_vertices']) *
                int(row['left_rooted_count']) * int(row['right_rooted_count']))
    assert expected == int(row['count']) == int(rows[k]['count'])

total = sum(int(r['count']) for r in rows.values())
audit = json.loads((cert / 'remote-audit.json').read_text())
assert total == int(audit['S1'])
lower = Fraction(total, 128*a2)
assert str(lower.numerator) == audit['strict_union_lower_numerator']
assert str(lower.denominator) == audit['strict_union_lower_denominator']

binary = out / 'ie_counter'
compiler = ['clang++', '-std=c++20', '-O3', '-Wall', '-Wextra', '-Werror',
            '-pedantic']
version = run(['clang++', '--version']).splitlines()[0]
run(compiler + [cert / 'curved_single_region_root.cpp', '-o', binary])
run([sys.executable, '-B', cert / 'gate_n5.py', binary,
     project / 'certificates/r042/evidence/diagonal.json', out / 'n5'])
for profile in ('3333333', '3434343', '4343434', '4444444'):
    target = out / (profile + '.json')
    run([binary, '7', target, profile])
    assert json.loads(target.read_text()) == rows[profile]
cross = json.loads((cert / 'subset-cross-audit.json').read_text())
assert all(r['match'] and r['cross_rooted_count'] == rows[r['profile']]['left_rooted_count']
           for r in cross['profiles'])
run([sys.executable, '-B', cert / 'binary_cut_transfer.py', out / 'transfer'])
result = dict(compiler=version, python=sys.version, profiles_checked=len(rows),
              n5_regenerated=32, n7_regenerated=4, n7_independent_stored=2,
              S1=total, A2=a2, union_lower=str(lower), children=records)
(out / 'audit.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
