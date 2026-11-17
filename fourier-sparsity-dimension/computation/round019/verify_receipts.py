#!/usr/bin/env python3
"""Audit R19 archived evidence, not a rerun of the exhaustive computations."""
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text())


def sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def monomial(indices):
    return sum(int(all((x >> i) & 1 for i in indices)) << x for x in range(64))


def main():
    m = read('macwilliams-rm36.json')
    assert m['source_sha256'] == sha('macwilliams_rm36.py')
    target = m['target_weight16_count']
    assert target == 232081500 and sum(m['target_weight_histogram']) == 2**42
    o = read('orbits-gpu.json')
    triples = list(itertools.combinations(range(6), 3))
    pairs = list(itertools.combinations(range(6), 2))
    assert o['cubic_basis_order'] == [''.join(f'x{i}' for i in t) for t in triples]
    assert o['quadratic_basis_order'] == [''.join(f'x{i}' for i in t) for t in pairs]
    assert [int(x, 16) for x in o['quadratic_truth_hex64']] == [monomial(t) for t in pairs]
    for row in o['orbit_representatives']:
        value = 0
        for i, t in enumerate(triples):
            if row['cubic_coefficient_mask'] >> i & 1:
                value ^= monomial(t)
        assert value == int(row['cubic_truth_hex64'], 16)
    # Independent Boolean substitution/Mobius check for the 60 generator basis images.
    for t in triples:
        for gen in range(3):
            values = []
            for x in range(64):
                bits = [(x >> j) & 1 for j in range(6)]
                mapped = ([bits[(j + 1) % 6] for j in range(6)] if gen == 0
                          else [bits[1], bits[0]] + bits[2:] if gen == 1
                          else [bits[0] ^ bits[1]] + bits[1:])
                values.append(int(all(mapped[j] for j in t)))
            for bit in range(6):
                for x in range(64):
                    if x >> bit & 1:
                        values[x] ^= values[x ^ (1 << bit)]
            actual = {tuple(j for j in range(6) if x >> j & 1)
                      for x in range(64) if values[x] and x.bit_count() == 3}
            if gen == 0:
                expected = {tuple(sorted((j + 1) % 6 for j in t))}
            elif gen == 1:
                expected = {tuple(sorted(1 if j == 0 else 0 if j == 1 else j for j in t))}
            else:
                expected = {t}
                if 0 in t and 1 not in t:
                    expected.add(tuple(sorted(1 if j == 0 else j for j in t)))
            assert actual == expected
    assert sum(r['orbit_size'] for r in o['orbit_representatives']) == 2**20
    c = read('full-result.json')
    r = read('full-run-receipt.json')
    assert r['hashes']['source_sha256'] == sha('audit_rm36_weight16.cpp')
    assert r['hashes']['runner_sha256'] == sha('run_rm36_full.py')
    assert r['hashes']['orbit_input_sha256'] == sha('orbits-gpu.json')
    assert r['result_sha256'] == sha('full-result.json')
    assert r['child_returncode'] == 0 and not r['timed_out'] and r['result'] == c
    assert c['total_weight16_count_weighted_by_orbit'] == target
    assert c['total_projectable_weight16_count_weighted_by_orbit'] == target
    assert c['total_nongraph_weight16_count_weighted_by_orbit'] == 0
    p = read('popcount-result.json')
    assert p['threads'] == 12 and p['weighted_weight16_count'] == target
    assert p['weighted_projectable_count'] == target and p['weighted_nongraph_count'] == 0
    assert sha('verify_rm36_popcount.cpp') in (HERE / 'LOCAL_CROSSCHECK_RECORD.md').read_text()
    assert sha('popcount-result.json') in (HERE / 'LOCAL_CROSSCHECK_RECORD.md').read_text()
    g = read('metal-result.json')
    assert g['input_sha256'] == sha('orbits-gpu.json')
    assert g['source_sha256'] == sha('R19.swift')
    assert g['status'] == 'passed' and g['weighted_weight16'] == target
    assert all(row['gpu_cpu_equal'] for row in g['rows'])
    for cpu, local, gpu, orbit in zip(c['orbits'], p['orbits'], g['rows'], o['orbit_representatives'], strict=True):
        assert cpu['representative'] == local['representative'] == gpu['mask'] == orbit['cubic_coefficient_mask']
        assert cpu['size'] == local['size'] == gpu['orbit_size'] == orbit['orbit_size']
        assert cpu['weight16_representative_count'] == local['weight16_count'] == gpu['fiber_weight16']
        assert cpu['projectable_representative_count'] == local['projectable_count'] == local['weight16_count']
    assert sum(row['weight16_representative_count'] for row in c['orbits']) == 4220
    d = read('d15-graph-projection-audit.json')
    assert d['source_sha256'] == sha('d15_graph_projection_audit.py')
    assert d['partitions_sha256'] == sha('../round006/weight8-partitions.json')
    assert d['audit_sha256'] == sha('../round006/common-mm-audit.json')
    assert (d['weight16_occurrences'], d['unique_weight16_switches'], d['nongraph_occurrences']) == (5680, 5076, 0)
    assert [row['graph_kernels'] for row in d['r18_positive_control']] == [16] * 4
    assert d['negative_control']['direct_graph_kernels'] == d['negative_control']['graph_kernels'] == 0
    n = read('nonstar-quotient-freeze-local-audit.json')
    assert n['source_sha256'] == sha('nonstar_quotient_freeze_local_audit.py')
    assert [row['j'] for row in n['rows']] == [2, 1, 1, 0]
    assert [row['local_affine_phase_image_rank'] for row in n['rows']] == [3, 2, 2, 1]
    print(json.dumps({'status': 'passed', 'scope': 'archived bindings, counts and generator basis action',
                      'weighted_weight16': target, 'representative_weight16': 4220,
                      'signed_cpu_gpu_coefficients': 12582912, 'generator_basis_images': 60}))


if __name__ == '__main__':
    main()
