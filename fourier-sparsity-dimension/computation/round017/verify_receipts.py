#!/usr/bin/env python3
"""Check archived source/input bindings and finite coverage, not general proofs."""
import ast
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / name).read_text())


def sha(name):
    return hashlib.sha256((ROOT / name).read_bytes()).hexdigest()


def main():
    for p in ROOT.glob('*.py'):
        ast.parse(p.read_text(), filename=p.name)
    for p in ROOT.glob('*.json'):
        json.loads(p.read_text())
    source = sha('verify_rank_two_families.py')
    n4, n5, n6 = [read(n) for n in ['local-n4.json', 'local-n5-anchored.json', 'local-n6-anchored.json']]
    for d in [n4, n5, n6]:
        assert d['source_sha256'] == source and d['failed_assertions'] == 0
    assert n4['affine_tested'] == 26387
    assert n4['subspace_level_counts'] == [1, 63, 651, 1395, 651, 63, 1]
    assert n5['subspace_level_counts'] == [1, 154, 616, 192]
    assert n6['subspace_level_counts'] == [1, 650, 3000, 1200, 48]
    remote = read('remote-enumeration-result.json')
    assert remote['source_sha256'] == sha('enumerate_affine_rank2_4x4.py')
    assert remote['status'] == 'passed'
    counts = [0] * 7
    for key, count in n4['valid_classifications'].items():
        counts[int(key.split('_')[0][1:])] += count
    assert counts == [35, 595, 840, 105, 0, 0, 0]
    assert counts == [remote['constant_rank_2_affine_families_by_dimension'][str(i)] for i in range(7)]
    assert sum(remote['enumerated_all_affine_subspaces_by_dimension'].values()) == 26387
    cases = read('nonparallel-cases.json')
    assert cases['source_sha256'] == sha('verify_nonparallel_blocks.py')
    assert cases['helper_sha256'] == source and cases['failed_assertions'] == 0
    assert [c['phase_count'] for c in cases['cases']] == [2048, 1024]
    assert sum(c['branch_checks'] for c in cases['cases']) == 13312
    gpu = read('nonparallel-metal-report.json')
    assert gpu['input_sha256'] == sha('nonparallel-cases.json')
    assert gpu['source_sha256'] == sha('NonparallelMetal.swift')
    assert gpu['status'] == 'passed' and gpu['total_phases'] == 3072
    assert gpu['all_command_buffers_completed'] and gpu['command_buffer_error'] == 'none'
    for a, b in zip(cases['cases'], gpu['cases']):
        assert a['name'] == b['name'] and a['dimension'] == b['dimension']
        assert a['phase_count'] == b['phases'] == b['cpu_phases']
        assert b['rank_ok'] and b['blocks_partition_ok'] and b['blocks_disjoint_ok']
        assert b['cpu_entrywise_equal'] and b['gpu_bad_abs8'] == 0
        offset = 0
        for c in b['chunks']:
            assert c['start'] == offset and c['bad_abs8'] == 0
            offset += c['count']
        assert offset == b['phases']
    manifest = dict((line.split()[1], line.split()[0]) for line in (ROOT / 'metal-manifest-sha256.txt').read_text().splitlines())
    for n in ['NonparallelMetal.swift', 'nonparallel-cases.json', 'nonparallel-metal-report.json']:
        assert manifest[n] == sha(n)
    cylinder = read('cylinder-audit.json')
    assert cylinder['schema'] == 'a28-r017-independent-cylinder-audit-v2'
    assert cylinder['input_sha256'] == sha('nonparallel-cases.json')
    assert cylinder['script_sha256'] == sha('check_cylinders.py')
    assert cylinder['failed_assertions'] == cylinder['exit_status'] == 0
    assert [c['freeze_rank_r'] for c in cylinder['cases']] == [2, 0]
    assert [c['V0_dimension'] for c in cylinder['cases']] == [9, 10]
    assert [c['normal_union_rank'] for c in cylinder['cases']] == [4, 5]
    for c in cylinder['cases']:
        assert c['disjoint_and_cover'] and c['walsh_support_slices_exact']
        assert c['cylinder_count'] == 8 and c['cylinder_codimension'] == 3
        assert c['parameter_projection_rank'] + c['K_dimension'] == c['V0_dimension']
        assert c['kernel_restriction_image_rank'] <= c['kernel_bound_m1_plus_2m2']
    print(json.dumps({'status': 'passed', 'unique_n4_families': sum(counts),
                      'local_n5_anchored_families': n5['affine_tested'],
                      'local_n6_anchored_families': n6['affine_tested'],
                      'full_CPU_GPU_phase_comparisons': gpu['total_phases'],
                      'scope': 'Source/input/result bindings and finite coverage; no general proof certification.'}))


if __name__ == '__main__':
    main()
