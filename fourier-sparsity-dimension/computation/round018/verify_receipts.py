#!/usr/bin/env python3
"""Check exact artifact bindings and finite coverage, not the general theorem."""
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
    data = read('graph-cases.json')
    assert data['source_sha256'] == sha('build_graph_cases.py')
    assert [x['dimension'] for x in data['cases']] == [9, 14]
    assert [x['phase_count'] for x in data['cases']] == [512, 16384]
    small = read('n4-double-switch-audit.json')
    assert small['source_sha256'] == sha('n4_double_switch_audit.py')
    assert small['total_bent_phases'] == 896 and small['weight_2q_switches'] == 390
    assert small['weight_2q_decomposable'] == 390 and small['weight_2q_nondecomposable'] == 0
    witness = read('nonaffine-four-block-audit.json')
    assert witness['source_sha256'] == sha('nonaffine_four_block_audit.py')
    assert witness['input_sha256'] == sha('graph-cases.json')
    assert witness['all_affine_q_flats'] == 11160
    assert witness['all_common_bent_phase_count'] == 512 and witness['all_common_bent_verified']
    assert witness['dual_supports_disjoint_partition']
    assert all(b['contained_all_affine_q_flats'] == 0 for b in witness['blocks'])
    assert int(witness['phase_mask_hex'], 16) == int(data['cases'][0]['f0_hex'], 16)
    cylinder = read('rank2-graph-cylinders.json')
    assert cylinder['source_sha256'] == sha('verify_rank2_graph_cylinders.py')
    assert cylinder['helper_sha256'] == sha('build_graph_cases.py')
    assert cylinder['joint_points_checked'] == 1024 and cylinder['failed_assertions'] == 0
    assert cylinder['rank2_compatible_with_fixed_Q_by_rank'] == {'2': 90, '4': 20, '6': 0}
    receipt = read('cpu-receipt.json')
    assert receipt['compile']['source_sha256'] == sha('verify_nonaffine_blocks.cpp')
    for r, n in enumerate(['cpu-r0-final.json', 'cpu-r1-final.json']):
        c, run = read(n), receipt['runs'][r]
        assert run['output_sha256'] == sha(n) and run['exit_code'] == 0
        assert c['phase_count_expected'] == c['phase_count_processed'] == data['cases'][r]['phase_count']
        assert c['all_checks_passed'] and c['generator_rank_ok']
        assert c['generator_rank'] == data['cases'][r]['dimension']
        for key in ['bad_phase_count', 'bad_walsh_term_count', 'maximum_abs_deviation',
                    'partial_support_failure_count', 'partial_support_size_failure_count',
                    'partial_nonzero_magnitude_failure_count', 'block_partition_failure_count',
                    'block_weight_failure_count']:
            assert c[key] == 0
        assert c['partial_block_checks'] == 2048
    gpu = read('graph-metal-report.json')
    assert gpu['input_sha256'] == sha('graph-cases.json')
    assert gpu['source_sha256'] == sha('GraphMetal.swift')
    assert gpu['status'] == 'passed' and gpu['phases'] == 16896
    assert len(gpu['rows']) == 2
    coefficients = 0
    for c, row in zip(data['cases'], gpu['rows']):
        assert row['n'] == c['n'] and row['q'] == c['q'] and row['phases'] == c['phase_count']
        assert row['rank_ok'] and row['cpu_equal'] and row['gpu_bad'] == 0
        coefficients += row['phases'] * (1 << row['n'])
    manifest = {line.split()[1]: line.split()[0] for line in (ROOT / 'metal-manifest-sha256.txt').read_text().splitlines()}
    for archived, remote in [('GraphMetal.swift', 'GraphMetal.swift'), ('graph-cases.json', 'graph-cases.json'), ('graph-metal-report.json', 'report.json')]:
        assert manifest[remote] == sha(archived)
    assert coefficients == 4227072
    print(json.dumps({'status': 'passed', 'nondecomposable_graph_blocks': 4,
                      'all_affine_3_flats_checked': 11160, 'full_CPU_GPU_phase_count': 16896,
                      'full_signed_coefficient_comparisons': coefficients,
                      'scope': 'Checks source/input/output bindings and finite coverage; no general mathematical proof certification.'}))


if __name__ == '__main__':
    main()
