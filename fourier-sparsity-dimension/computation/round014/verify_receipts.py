#!/usr/bin/env python3
"""Check archived bindings and certificates, without repeating the remote scans.

This validates receipts; it does not re-execute the GPU or the 29,085 DPs.
"""
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / name).read_text())


def sha(name):
    return hashlib.sha256((ROOT / name).read_bytes()).hexdigest()


def main():
    for path in ROOT.glob('*.json'):
        json.loads(path.read_text())
    for path in ROOT.glob('*.py'):
        ast.parse(path.read_text(), filename=str(path))

    gen = read('remote-generation-summary.json')
    cpu = read('remote-cpu-summary.json')
    audit = read('remote-parent-audit.json')
    assert gen['source_sha256'] == sha('build_addback_tasks.py')
    assert gen['base_codes'] == 2825
    assert gen['all_algebra_pairs'] == 52897
    assert gen['stable_pairs_proved_without_search'] == 23812
    assert gen['tasks'] == cpu['successful'] == cpu['records'] == cpu['tasks_total'] == audit['coverage'] == 29085
    assert gen['tasks_sha256'] == cpu['tasks_sha256'] == audit['tasks_sha256']
    assert audit['summary_sha256'] == sha('remote-cpu-summary.json')
    assert audit['source_sha256'] == sha('verify_batch_results.py')
    assert audit['results_sha256'] == cpu['results_sha256']
    assert cpu['unstarted_tasks'] == 0 and not cpu['failures'] and not cpu['counterexamples']
    boundary = read('boundary-partition-audit.json')
    assert boundary['input_sha256'] == sha('remote-cpu-summary.json')
    assert len(boundary['records']) == 20
    assert all(sum(r['partition_terms']) == r['kappa'] == r['upper'] for r in boundary['records'])

    inputs = read('selector-input.json')['records']
    gpu = read('metal-walsh-report.json')
    results = [json.loads(s) for s in (ROOT / 'selector-capacity-results.jsonl').read_text().splitlines()]
    by_id = {r['id']: r for r in results}
    assert len(inputs) == len(gpu['rows']) == len(by_id) == len(results) == 115
    assert gpu['source_sha256'] == sha('MetalWalshVerifier.swift')
    assert gpu['input_sha256'] == sha('selector-input.json')
    assert gpu['backend'] == 'Swift Metal int32 threadgroup FWHT'
    assert gpu['all_command_buffers_completed'] and gpu['total_bad_abs8'] == 0
    assert gpu['calibration']['equal']
    for i, (inp, row) in enumerate(zip(inputs, gpu['rows'])):
        r = by_id[i]
        assert row['index'] == i and row['rank'] == inp['rank'] == r['D']
        assert r['task']['lift_basis'] == inp['basis_masks']
        assert row['phases'] == 1 << inp['rank']
        assert row['bad_abs8'] == 0 and row['cpu_equal'] and row['stabilizer_count'] == 2
        assert r['returncode'] == 0 and r['status'] == 'ok' and r['kappa'] <= 16
    assert sum(r['phases'] for r in gpu['rows']) == 1347584
    assert sum(r['cpu_phases'] for r in gpu['rows']) == gpu['total_cpu_phases'] == 33680
    assert Counter(r['kappa'] for r in results) == {10: 8, 11: 24, 16: 83}
    selector_audit = read('selector-result-audit.json')
    assert selector_audit['results_sha256'] == sha('selector-capacity-results.jsonl')
    assert selector_audit['summary_sha256'] == sha('selector-capacity-summary.json')
    bind = read('cpu-gpu-binding-audit.json')
    for key, name in [('GPU_source_sha256', 'MetalWalshVerifier.swift'),
                      ('GPU_report_sha256', 'metal-walsh-report.json'),
                      ('input_sha256', 'selector-input.json'),
                      ('CPU_results_sha256', 'selector-capacity-results.jsonl')]:
        assert bind[key] == sha(name)

    remote = read('remote-eight-point-summary.json')
    assert remote['source_sha256'] == sha('capacity_dp_parallel.cpp')
    assert remote['input_sha256'] == sha('eight-point-candidates.json')
    assert remote['driver_sha256'] == sha('run_d18_remote.py')
    assert remote['validation_passed'] and remote['all_authorized_runs_ok']
    assert len(remote['runs']) == 4
    assert all(r['returncode'] == 0 and not r['timed_out'] for r in remote['runs'])
    for stem in ['local-eight-point', 'remote-eight-point-free-pair', 'remote-eight-point-mixed']:
        receipt = read(stem + '-verification.json')
        result = read(stem + '-result.json')
        assert receipt['source_sha256'] == sha('verify_eight_point_case.py')
        assert receipt['input_sha256'] == sha('eight-point-candidates.json')
        assert receipt['result_sha256'] == sha(stem + '-result.json')
        assert receipt['result'] == result
        assert sum(receipt['partition_terms']) == result['kappa'] == 24
        assert sum(c['kappa'] for c in receipt['direct_sum_components']) == 24
        assert all(c['serial_exit'] == c['parallel_exit'] == 0 for c in receipt['direct_sum_components'])
    for r, stem in zip(remote['runs'][2:], ['remote-eight-point-free-pair', 'remote-eight-point-mixed']):
        # The compact remote certificates retain only the mathematical output.
        compact = read(stem + '-result.json')
        assert {k: r['parsed'][k] for k in compact} == compact
    mm = read('mm-affine-jump.json')
    assert mm['status'] == 'passed' and mm['source_sha256'] == sha('verify_mm_affine_capacity_jump.py')
    print(json.dumps({'status': 'passed', 'scope': 'Archived input/source/result bindings and certificate summaries; not a fresh GPU or full DP run.', 'remote_tasks': 29085, 'selector_codes': 115, 'eight_point_codes': 3}))


if __name__ == '__main__':
    main()
