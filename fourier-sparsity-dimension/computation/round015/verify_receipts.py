#!/usr/bin/env python3
"""Check final source bindings, independently reconstruct Q4 bases and validate receipts."""
import ast
import hashlib
import json
from pathlib import Path
from triangle_graph_scan import basis, words

HERE = Path(__file__).resolve().parent


def sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def read(name):
    return json.loads((HERE / name).read_text())


def main():
    for f in HERE.glob('*.json'):
        json.loads(f.read_text())
    for f in HERE.glob('*.py'):
        ast.parse(f.read_text(), filename=str(f))
    old = read('../round014/remote-parent-audit.json')
    remote = read('remote-core-capacity-summary.json')
    assert remote['source_sha256'] == sha('audit_core_capacity.py')
    assert remote['tasks_sha256'] == old['tasks_sha256']
    assert remote['results_sha256'] == old['results_sha256']
    assert remote['status'] == 'passed' and remote['tasks'] == 29085
    assert remote['gap_histogram'] == {'0': 26565, '1': 1680, '2': 840}
    assert remote['violations'] == 0 and remote['strict'] == 2520
    assert remote['counterexample_sha256'] == sha('remote-counterexample-result.json')
    strict = read('strict-examples.json')
    assert strict['source_sha256'] == sha('extract_strict_examples.py')
    for k in ['comparison_sha256', 'tasks_sha256', 'results_sha256']:
        assert strict['sources'][k] == remote[k]
    checked = read('strict-examples-verification.json')
    assert checked['source_sha256'] == sha('verify_strict_examples.py')
    assert checked['utility_sha256'] == sha('triangle_graph_scan.py')
    assert checked['input_sha256'] == sha('strict-examples.json')
    assert [r['gap'] for r in checked['rows']] == [1, 2]
    for r in checked['rows']:
        ex = strict['examples'][r['name']]
        assert r['exit_code'] == 0 and r['T']['kappa'] == ex['kappa_T']
        assert r['C']['kappa'] == ex['kappa_C'] and r['B']['kappa'] == ex['kappa_B']
    kernel = read('nonunital-kernel-verification.json')
    assert kernel['source_sha256'] == sha('verify_nonunital_kernel.py')
    assert kernel['status'] == 'passed' and kernel['cases'] == 16
    assert all(r['sigma_N'] == r['kappa_B'] for r in kernel['rows'])
    tri, res = read('triangle-scan.json'), read('triangle-resource.json')
    assert tri['source_sha256'] == res['source_sha256'] == sha('triangle_graph_scan.py')
    assert tri['checked'] == 16384 and tri['unstable_E_cases'] == 256
    assert sum(tri['gap_histogram'].values()) == tri['checked']
    assert max(map(int, tri['gap_histogram'])) < 0 and res['exit_code'] == 0

    inp, gpu = read('q4-input.json'), read('q4-metal-report.json')
    assert gpu['source_sha256'] == sha('MetalWalshVerifier_r015.swift')
    assert gpu['input_sha256'] == sha('q4-input.json')
    assert gpu['status'] == 'passed' and gpu['all_command_buffers_completed']
    assert gpu['total_bad_abs8'] == 0 and gpu['total_cpu_phases'] == 98304
    manifest = {line.split()[1]: line.split()[0]
                for line in (HERE / 'q4-metal-manifest-sha256.txt').read_text().splitlines()}
    for name in ['MetalWalshVerifier_r015.swift', 'q4-input.json', 'q4-metal-report.json']:
        assert manifest[name] == sha(name)
    assert len(inp['records']) == len(gpu['rows']) == 3
    full = (1 << 64)-1
    coords = [sum(((i >> j) & 1) << i for i in range(64)) for j in range(6)]
    x, z0, y, w0, z1, w1 = coords
    for row, result, a in zip(inp['records'], gpu['rows'], [x, y, x ^ y]):
        # Reconstruct using Boolean truth masks and bitwise products,
        # independently of the generator's pointwise lambda construction.
        generators = []
        for u in [w0, full ^ w0]:
            for v in [w1, full ^ w1]:
                generators += [u & v, u & v & a]
        generators += [z0, z0 & a]
        for u in [w0, full ^ w0]:
            generators += [z1 & u, z1 & u & a]
        sb = basis(generators)
        tb = basis([*sb, full, x, y])
        assert len(sb) == row['S_rank'] == 14
        assert len(tb) == row['rank'] == result['rank'] == 15
        assert tb == basis(row['basis_masks'])
        assert all(len(basis([*tb, g])) == 15 for g in [full, *coords])
        assert result['phases'] == result['cpu_phases'] == 32768
        assert result['bad_abs8'] == 0 and result['cpu_equal']
        assert result['stabilizer_count'] == 2
    print(json.dumps({'status':'passed','scope':'Receipt/source bindings and independent Q4 basis reconstruction; does not rerun the historical T DPs or Metal kernel.', 'Q4_ranks':[15,15,15], 'full_CPU_GPU_compared_phases':98304}))


if __name__ == '__main__':
    main()
