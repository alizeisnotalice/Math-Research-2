#!/usr/bin/env python3
"""Compare all h/E/far residues and CRT values to independent joint tables."""
import argparse
import hashlib
import json
from math import factorial, prod
from pathlib import Path

PRIMES = [9223372036854775783, 9223372036854775643,
          9223372036854775549, 9223372036854775507]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--runs', type=Path, required=True)
    p.add_argument('--references', type=Path, required=True)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    manifest = json.loads(a.manifest.read_text())
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    assert sha(a.source) == manifest['source_sha256']
    modulus = prod(PRIMES)
    checks = []
    for ns, vertices in manifest['gate_points'].items():
        n = int(ns)
        assert modulus > 2*(n-1)**2*factorial(n*n-1)
        for v in vertices:
            name = 'gap-n3-all.json' if n == 3 else f'gap-n5-v{v}.json'
            ref = a.references/name
            assert sha(ref) == manifest['reference_sha256'][name]
            document = json.loads(ref.read_text())
            assert document['n'] == n
            hole = next(x for x in document['holes'] if x['v'] == v)
            expected = {}
            for t in range(1, 6):
                cells = [(d, e, int(c)) for L, d, e, c in hole['joint_H_L_d_e2'] if L == t]
                expected[f'h{t}_residue'] = sum(c for d, e, c in cells)
                assert expected[f'h{t}_residue'] == int(hole['histogram_L'][t-1])
                expected[f'E{t}_residue'] = sum(e*c for d, e, c in cells)
                for k in range(1, 4):
                    expected[f'far{k}_h{t}_residue'] = sum(c for d, e, c in cells if 2*d > k*n)
            rows = []
            for pi, prime in enumerate(PRIMES):
                raw = a.runs/f'n{n}-v{v}-p{pi}.json'
                row = json.loads(raw.read_text())
                assert (row['n'], row['v'], row['prime_index'], int(row['prime'])) == (n, v, pi, prime)
                assert row['algorithm'] == 'gap-spatial-h1-to-h5-modular-layer-rolling'
                assert row['max_modeled_payload_bytes'] <= row['payload_limit_bytes']
                for key, count in expected.items():
                    assert 0 <= int(row[key]) < prime
                    assert int(row[key]) == count % prime, (n, v, pi, key, row[key], count)
                rows.append(row)
            for key, count in expected.items():
                direct = sum(int(row[key])*(modulus//prime)*pow(modulus//prime, -1, prime)
                             for row, prime in zip(rows, PRIMES)) % modulus
                value, factor = 0, 1
                for row, prime in zip(rows, PRIMES):
                    value += ((int(row[key])-value)*pow(factor, -1, prime) % prime)*factor
                    factor *= prime
                assert direct == value == count, (n, v, key)
            checks.append({'n': n, 'v': v, 'exact_fields': expected})
    assert len(checks) == 15
    result = {'status': 'passed', 'points': 15, 'residue_files': 60,
              'integer_fields_per_point': 25, 'checks': checks,
              'source_sha256': sha(a.source), 'manifest_sha256': sha(a.manifest),
              'scope': 'n3 all holes and n5 six D4 classes: exact joint-table residue and dual-CRT comparison; no n7 production claim'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['status', 'points', 'residue_files', 'integer_fields_per_point']}))


if __name__ == '__main__':
    main()
