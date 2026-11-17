#!/usr/bin/env python3
"""CPU-only rejection tests for the d15 resume adapter's state primitives."""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
from pathlib import Path

import torch


HERE = Path(__file__).resolve().parent
ADAPTER = HERE / 'resume_gpu_factor_search_d15.py'


def reject(label, fn):
    try:
        fn()
    except (ValueError, SystemExit):
        return label
    raise AssertionError(f'accepted invalid {label}')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--state', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('d15_resume', ADAPTER)
    assert spec and spec.loader
    adapter = importlib.util.module_from_spec(spec); spec.loader.exec_module(adapter)
    state = torch.load(args.state, map_location='cpu', weights_only=True)
    m, n, d, starts = 64, 104, 15, 49152
    free_shape, theta_shape = (starts, m + n - 2*d - 2), (starts, 1)
    expected = (state['step'] % 2000) + 1
    adapter.finite_cpu_tensor(state['free'], free_shape, 'actual-free')
    adapter.validate_optimizer(state['optimizer'], free_shape, theta_shape, expected)
    bad_shape = state['free'][:1]
    bad_opt = {'state': state['optimizer']['state'], 'param_groups': [copy.copy(state['optimizer']['param_groups'][0])]}
    bad_opt['param_groups'][0]['lr'] = .1
    out = {
        'status': 'passed_cpu_rejection_tests', 'state_step': state['step'],
        'rejected': [
            reject('free_wrong_shape', lambda: adapter.finite_cpu_tensor(bad_shape, free_shape, 'bad-free')),
            reject('adam_wrong_lr', lambda: adapter.validate_optimizer(bad_opt, free_shape, theta_shape, expected)),
        ],
        'claim_boundary': 'Tests loader primitive rejection only; no MPS run or resume was started.',
    }
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps(out, sort_keys=True))


if __name__ == '__main__':
    main()
