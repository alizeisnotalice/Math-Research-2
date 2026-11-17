#!/usr/bin/env python3
"""Verify input/source/result bindings and complete CPU/GPU family coverage."""
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
    for f in ROOT.glob('*.json'):
        json.loads(f.read_text())
    for f in ROOT.glob('*.py'):
        ast.parse(f.read_text(), filename=str(f))
    seed, checked = read('extension-seed.json'), read('seed-verification.json')
    assert checked['source_sha256'] == sha('build_extension_seed.py')
    assert checked['seed_sha256'] == sha('extension-seed.json')
    assert checked['bent_truth_tables'] == 896
    assert checked['valid_quotient_subspaces_by_dimension'] == [1, 27, 45]
    assert checked['inclusion_maximal_base_extensions_by_dimension'] == {'5': 0, '6': 0, '7': 45}
    assert checked['maximum_base_extension_dimension'] == 7
    assert checked['maximum_base_extension_capacities'] == [8]*45
    assert checked['selector_small_test']['exactly_selector']
    assert checked['selector_small_test']['common_legal_directions'] == 128
    assert checked['constant_only_rank2_reduction']['branch_spectra_checked'] == 256
    cap = read('base-capacity-audit.json')
    assert cap['source_sha256'] == sha('verify_base_capacities.py')
    assert cap['utility_sha256'] == sha('build_extension_seed.py')
    assert cap['input_sha256'] == sha('seed-verification.json')
    assert cap['cases'] == len(cap['rows']) == 45
    assert [r['index'] for r in cap['rows']] == list(range(45))
    assert all(r['kappa'] == 8 and r['unital_subspaces'] == 2825 for r in cap['rows'])
    assert cap['total_unital_subspaces'] == 127125 and cap['total_algebras'] == 4335

    masks = list(map(int, (ROOT / 'seed_masks.txt').read_text().split()))
    assert masks == seed['relative_bent_phase_masks']
    assert len(masks) == len(set(masks)) == seed['seed_count'] == 896
    assert all(0 <= v < 65536 for v in masks) and 0 in masks
    extracted = read('seed-extraction.json')
    assert extracted['extractor_sha256'] == sha('extract_seed_masks.py')
    assert extracted['json_sha256'] == sha('extension-seed.json')
    assert extracted['text_sha256'] == sha('seed_masks.txt')
    total = 2*len(masks)**2
    assert seed['product_family_count'] == total == 1605632
    cpu, receipt = read('cpu-extension-result.json'), read('cpu-extension-receipt.json')
    assert receipt['child_returncode'] == 0 and not receipt['timed_out']
    assert receipt['result'] == cpu and receipt['result_sha256'] == sha('cpu-extension-result.json')
    for name, key in [('verify_extension_cpu.cpp','cpp_sha256'), ('run_extension_cpu.py','runner_sha256'),
                      ('extension-seed.json','seed_json_sha256'), ('seed_masks.txt','seed_text_sha256'),
                      ('seed-extraction.json','seed_extraction_sha256')]:
        assert receipt['hashes'][key] == sha(name)
    assert cpu['phase_count_expected'] == cpu['phase_count_processed'] == total
    assert cpu['full_coverage'] and cpu['all_abs_walsh_equal_8']
    assert cpu['bad_phase_count'] == cpu['bad_walsh_term_count'] == cpu['maximum_bad_abs_deviation_from_8'] == 0

    gpu = read('extension-metal-report.json')
    assert gpu['source_sha256'] == sha('ExtensionMetal.swift')
    assert gpu['input_sha256'] == sha('extension-seed.json')
    assert gpu['status'] == 'passed' and gpu['all_command_buffers_completed']
    assert gpu['command_buffer_error'] == 'none' and gpu['bad_abs8'] == 0
    assert gpu['records_expected'] == gpu['records_completed'] == total
    assert gpu['chunks_expected'] == gpu['chunks_completed'] == len(gpu['chunks']) == 784
    assert gpu['chunk_phases'] == gpu['calibration']['phases'] == 2048
    assert gpu['calibration']['equal']
    start = 0
    for i, r in enumerate(gpu['chunks']):
        assert r['index'] == i and r['start'] == start and r['count'] == 2048
        assert r['bad_abs8'] == 0
        start += r['count']
    assert start == total
    manifest = {s.split()[1]: s.split()[0] for s in (ROOT / 'metal-manifest-sha256.txt').read_text().splitlines()}
    for name in ['ExtensionMetal.swift', 'extension-seed.json', 'extension-metal-report.json']:
        assert manifest[name] == sha(name)
    print(json.dumps({'status':'passed', 'seed_completions':45,
                      'six_variable_maximal_completions_from_theorem':45**2,
                      'CPU_full_phases':total, 'GPU_full_phases':total,
                      'CPU_GPU_coefficientwise_calibration_phases':2048,
                      'scope':'Checks source/input/output bindings and coverage; full-family CPU/GPU coefficientwise equality is not asserted. Completeness among arbitrary extra directions uses proof016.'}))


if __name__ == '__main__':
    main()
