#!/usr/bin/env python3
"""Make a dimension-checked Metal input for exactly one GPU52 task.

This is deliberately independent of the task-2 prototype: every shape used by
the Swift program is copied from the selected task and checked before writing.
"""
import argparse
import hashlib
import json
from pathlib import Path

SOURCE = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-macbook-gpu52-direct/inputs.json')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_payload(task_number):
    doc = json.loads(SOURCE.read_text())
    tasks, dc, ds = doc['tasks'], doc['DC'], doc['DS']
    if not (0 <= task_number < len(tasks)):
        raise ValueError(f'task {task_number} outside 0..{len(tasks)-1}')
    task = tasks[task_number]
    if task.get('id') != task_number:
        raise ValueError(f'task index/id disagreement: {task.get("id")} != {task_number}')
    nfreq, ncoef = len(dc), len(dc[0]) if dc else 0
    if not (1 <= nfreq <= 64 and 1 <= ncoef <= 53):
        raise ValueError(f'unsupported Metal dimensions {(nfreq, ncoef)}')
    if any(len(r) != ncoef for r in dc) or len(ds) != nfreq or any(len(r) != ncoef for r in ds):
        raise ValueError('DC/DS are not rectangular with identical dimensions')
    groups = task['groups']
    if not (1 < groups < 32):
        raise ValueError(f'groups={groups} does not fit a UInt32 non-control mask')
    if len(task['LR']) != groups or len(task['LI']) != groups:
        raise ValueError('LR/LI group dimension mismatch')
    if any(len(r) != nfreq for r in task['LR']) or any(len(r) != nfreq for r in task['LI']):
        raise ValueError('LR/LI frequency dimension mismatch')
    masks = int(task['masks'])
    if not (0 < masks <= (1 << (groups - 1))):
        raise ValueError(f'bad task mask count {masks}')
    return {
        'task': task_number, 'groups': groups, 'masks': masks,
        'LR': task['LR'], 'LI': task['LI'], 'DC': dc, 'DS': ds,
        'source_inputs_sha256': sha(SOURCE),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--task', type=int, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    data = checked_payload(a.task)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, 'w') as f:
        f.write(json.dumps(data, sort_keys=True, separators=(',', ':')) + '\n')
        f.flush()
        import os
        os.fsync(f.fileno())
    print(json.dumps({'task': a.task, 'out': str(a.out), 'sha256': sha(a.out),
                      'groups': data['groups'], 'nfreq': len(data['DC']),
                      'ncoef': len(data['DC'][0]), 'masks': data['masks']}, sort_keys=True))


if __name__ == '__main__':
    main()
