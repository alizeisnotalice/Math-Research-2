#!/usr/bin/env python3
"""Independent exact terminal-tree check of the R24 affine family certificate.

Uses Fraction and JSON only: no family kernel is imported and no halo geometry
is recomputed. The archived audit records supply earlier independent tiny/full-
label checks. Source identities are checked here against actual files.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(snapshot, contract, baseline, parallel, family, reference, certificate,
          envelope, driver, halo, proof, helper):
    raw = snapshot.read_bytes()
    record = json.loads(raw)
    c = json.loads(contract.read_text())
    b = json.loads(baseline.read_text())
    p = json.loads(parallel.read_text())
    ref = json.loads(reference.read_text())
    cert = json.loads(certificate.read_text())
    source_hashes = {str(x): sha(x) for x in
                     [contract, baseline, parallel, family, reference,
                      certificate, envelope, driver, halo, proof, helper]}
    # Link saved final data to prior independent checks and actual source bytes.
    assert record['input_sha256'] == c['input_sha256'] == sha(family)
    assert record['driver_sha256'] == p['source_sha256'] == sha(driver)
    for suffix, actual in [('mass_family_envelope.py', envelope),
                           ('periodic_halo.py', halo),
                           ('affine-mass-rectangle-enclosure.md', proof),
                           ('check_affine_mass_enclosure.py', helper)]:
        expected = [v for k, v in b['sources'].items() if k.endswith('/' + suffix)]
        assert expected == [sha(actual)]
    assert record['envelope_sha256'] == sha(envelope)
    assert p['baseline_sha256'] == sha(baseline)
    assert record['reference_input_sha256'] == sha(reference) == cert['input_sha256']
    assert record['reference_certificate_sha256'] == sha(certificate)
    assert cert['code_sha256'] == sha(halo)
    assert F(record['lower']) == F(cert['inner_ratio'])
    assert F(cert['mass']) == F(c['actual_constant_mass'])
    assert F(cert['inner_ratio']) == F(cert['inner_area']) * F(ref['lam']) / F(cert['mass'])

    data = json.loads(family.read_text())
    a0, a1 = map(F, c['bounds']['a'])
    b0, b1 = map(F, c['bounds']['b'])
    assert list(map(F, data['a_bounds'])) == [a0, a1]
    assert list(map(F, data['b_bounds'])) == [b0, b1]
    a, z = map(F, record['reference_point'])
    assert a0 <= a <= a1 and b0 <= z <= b1
    bx, by, lam = (F(data['base'][k]) for k in ('Px', 'Py', 'lam'))
    nx, ny = data['copy_factors']
    assert (F(ref['Px']), F(ref['Py']), F(ref['lam'])) == (nx * bx, ny * by, lam)
    assert len(ref['pts']) == len(c['all_label_records']) == 204
    for row, label in zip(ref['pts'], c['all_label_records']):
        x, y, w = map(F, row)
        assert [str(x), str(y)] == label['coordinate']
        assert w == F(label['intercept']) + a * F(label['u']) + z * F(label['v'])
    assert record['family']['all_corner_labels_strictly_positive']
    assert c['all_corners_strictly_positive'] and F(c['minimum_corner_mass']) > 0
    assert c['sum_u'] == c['sum_v'] == '0'
    for key, contract_key in [
        ('mass', 'actual_constant_mass'),
        ('minimum_corner_label_mass', 'minimum_corner_mass'),
        ('uniform_kappa', 'uniform_kappa'), ('delta', 'delta'),
        ('actual_maximum_distance', 'maximum_lift_distance'),
        ('radius', 'R'), ('uniform_tail_bound', 'uniform_tail'),
        ('tail_margin', 'margin')]:
        assert F(record['family'][key]) == F(c[contract_key])
    assert F(record['family']['uniform_tail_bound']) < lam
    bits = record['family']['bits']
    assert type(bits) is int and bits >= 0
    assert F(record['family']['dyadic_outward_radius_error_at_most']) == F(1, 1 << bits)

    nodes = record['nodes']
    assert nodes and nodes[0]['parent'] is None
    assert tuple(map(F, nodes[0]['cell'])) == (a0, a1, b0, b1)
    children = []
    splits = 0
    for i, node in enumerate(nodes):
        assert type(node['id']) is int and node['id'] == i
        al, ah, bl, bh = map(F, node['cell'])
        assert a0 <= al < ah <= a1 and b0 <= bl < bh <= b1
        assert node['status'] in ('split', 'active', 'pruned')
        assert F(node['upper']) >= 0
        if i:
            assert type(node['parent']) is int and 0 <= node['parent'] < i
            assert F(node['upper']) <= F(nodes[node['parent']]['upper'])
        if node['status'] == 'split':
            splits += 1
            ids = node['children']
            assert len(ids) == 2 and len(set(ids)) == 2
            assert all(type(j) is int and i < j < len(nodes) for j in ids)
            left, right = [nodes[j] for j in ids]
            assert left['parent'] == right['parent'] == i
            lc, rc = tuple(map(F, left['cell'])), tuple(map(F, right['cell']))
            am, bm = (al + ah) / 2, (bl + bh) / 2
            assert (lc, rc) in [((al, am, bl, bh), (am, ah, bl, bh)),
                                ((al, ah, bl, bm), (al, ah, bm, bh))]
            children.extend(ids)
        else:
            assert node['children'] == []
    assert sorted(children) == list(range(1, len(nodes)))
    assert len(nodes) == 2 * splits + 1
    leaves = [n for n in nodes if n['status'] != 'split']
    assert len(leaves) == splits + 1 == record['leaf_count']
    assert splits == record['splits']
    assert len(leaves) <= record['max_leaves']
    domain_area = (a1 - a0) * (b1 - b0)
    total_area = sum(((F(n['cell'][1]) - F(n['cell'][0])) *
                      (F(n['cell'][3]) - F(n['cell'][2])) for n in leaves), F(0))
    assert total_area == domain_area
    # Every node has exactly one earlier parent, all nonroot nodes occur as
    # children exactly once, and every split is an exact adjacent midpoint
    # split. Induction proves full coverage and disjoint leaf interiors.
    lower, upper, gap, target = map(F, [record[k] for k in
                                      ('lower', 'upper', 'gap', 'target_gap')])
    assert lower >= 0 and target > 0
    assert all(F(n['upper']) < lower for n in leaves if n['status'] == 'pruned')
    assert all(F(n['upper']) >= lower for n in leaves if n['status'] == 'active')
    maximum = max([lower] + [F(n['upper']) for n in leaves])
    assert maximum == upper and upper - lower == gap
    assert record['status'] == 'leaf_limit'
    assert len(leaves) == record['max_leaves'] == 5000 and gap > target
    maximizing_leaves = [n['id'] for n in leaves if F(n['upper']) == upper]
    assert maximizing_leaves
    return {
        'status': 'passed',
        'scope': 'terminal exact tree, full-label reference binding and source identity; no real-family per-box or per-leaf geometry replay',
        'snapshot_sha256': hashlib.sha256(raw).hexdigest(),
        'snapshot_bytes': len(raw), 'snapshot_status': record['status'],
        'nodes': len(nodes), 'splits': splits, 'leaves': len(leaves),
        'leaf_status_counts': dict(Counter(n['status'] for n in leaves)),
        'all_binary_splits_and_backlinks_checked': True,
        'full_domain_cover_and_interior_disjointness_proved_by_tree': True,
        'sum_leaf_areas': str(total_area),
        'reference_point': record['reference_point'],
        'reference_all_labels_exactly_in_family': True,
        'reference_certificate_identity_and_mass_ratio_checked': True,
        'family_bounds': c['bounds'], 'labels': 204,
        'lower': str(lower), 'upper': str(upper), 'gap': str(gap),
        'target_gap': str(target), 'target_gap_reached': gap <= target,
        'maximum_upper_leaf_ids': maximizing_leaves,
        'bits': bits, 'source_hashes': source_hashes,
        'not_recomputed': ['real 5081629-box enumeration',
                           'real 24646-box retained list',
                           'real per-node rectangle unions',
                           'R23 reference halo area']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['snapshot', 'contract', 'baseline', 'parallel', 'family',
                'reference', 'certificate', 'envelope', 'driver', 'halo',
                'proof', 'helper', 'output']:
        parser.add_argument('--' + key, required=True, type=Path)
    args = {key: value.resolve() for key, value in vars(parser.parse_args()).items()}
    output = args.pop('output')
    start = time.process_time()
    result = check(**args)
    result.update(python=sys.version, platform=platform.platform(),
                  cpu_seconds=time.process_time() - start,
                  script_sha256=sha(Path(__file__)),
                  actual_command=[sys.executable, str(Path(__file__).resolve())] +
                  [item for key, value in [*args.items(), ('output', output)]
                   for item in ('--' + key, str(value))])
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('source_hashes', 'actual_command')}, indent=2))


if __name__ == '__main__':
    main()
