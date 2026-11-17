#!/usr/bin/env python3
"""Extract canonical strict R15 examples from completed, hash-bound R14 data."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tasks', required=True)
    ap.add_argument('--results', required=True)
    ap.add_argument('--comparison', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    tasks = {x['id']: x for x in rows(args.tasks)}
    results = {x['id']: x for x in rows(args.results)}
    comparison = rows(args.comparison)
    if set(tasks) != set(results) or set(tasks) != {x['id'] for x in comparison}:
        raise ValueError('ID coverage mismatch')
    out = {'schema': 'a28-r015-strict-example-extract-v1',
           'ordering': '(lift_dimension,id)',
           'sources': {'tasks_sha256': sha(args.tasks), 'results_sha256': sha(args.results),
                       'comparison_sha256': sha(args.comparison)}, 'examples': {}}
    for gap in (1, 2):
        candidates = [x for x in comparison if x['gap'] == gap]
        if not candidates:
            raise ValueError(f'no gap={gap} row')
        chosen = min(candidates, key=lambda x: (x['lift_dimension'], x['id']))
        i = chosen['id']; task = tasks[i]; result = results[i]
        if result['task'] != task or result['kappa'] != chosen['kappa_t']:
            raise ValueError('task/result binding mismatch')
        out['examples'][f'gap_{gap}'] = {
            'comparison': chosen,
            'C_basis': task['C_basis'], 'B_basis': task['B_basis'], 'E_cells': task['E_cells'],
            'T_lift_basis': task['lift_basis'], 'T_lift_dimension': task['lift_dimension'],
            'kappa_C': task['base_capacity'], 'kappa_B': chosen['kappa_b'],
            'kappa_T': result['kappa'], 'rhs_kappa_C_plus_kappa_B': chosen['rhs'],
            'T_optimal_partition': result['partition'],
            'completed_T_returncode': result['returncode'],
            'completed_T_status': result['status'],
        }
    out['source_sha256'] = sha(__file__)
    Path(args.out).write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'gap_1_id': out['examples']['gap_1']['comparison']['id'],
                      'gap_2_id': out['examples']['gap_2']['comparison']['id']}))


if __name__ == '__main__':
    main()
