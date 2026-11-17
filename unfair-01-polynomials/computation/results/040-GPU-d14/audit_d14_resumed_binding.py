#!/usr/bin/env python3
"""Bind a completed d14 resume to its immutable inputs and terminal records.

This complements (rather than replaces) verify_prefix_fraction.py: that
program replays the exact saved coefficient rows, while this audit checks that
the replay was applied to the declared resumed job and source inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import torch


ORIGINAL_SHA = 'f0044bf5a0667b1825d1b58662243bc6c9b49682c809d045680ef5c53d348d9f'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok: bool, message: str) -> None:
    if not ok:
        raise SystemExit(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--resume-worker', type=Path, required=True)
    parser.add_argument('--original-worker', type=Path, required=True)
    parser.add_argument('--fraction-audit', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run = args.run.resolve()
    config = json.loads((run / 'config.json').read_text())
    provenance = json.loads((run / 'resume-provenance.json').read_text())
    release = json.loads((run / 'release.json').read_text())
    checkpoint = json.loads((run / 'checkpoint.json').read_text())
    fraction = json.loads(args.fraction_audit.read_text())
    require(fraction['status'] == 'saved_prefix_samples_fraction_and_pool_audit_passed', 'Fraction audit did not pass')
    require(config['script_sha256'] == sha(args.resume_worker), 'config does not bind audited resume worker')
    require(config['original_worker_sha256'] == ORIGINAL_SHA == sha(args.original_worker), 'original worker source binding mismatch')
    require(provenance['resume_worker_sha256'] == sha(args.resume_worker), 'provenance resume worker mismatch')
    require(provenance['source_worker_sha256'] == ORIGINAL_SHA, 'provenance original worker mismatch')
    old_run = Path(config['resume_run'])
    old_state = Path(config['resume_state'])
    require(old_run.is_dir() and old_state.is_file(), 'declared resume input no longer readable')
    require(provenance['resume_run'] == str(old_run.resolve()), 'resume run binding mismatch')
    require(provenance['resume_state'] == str(old_state.resolve()), 'resume state binding mismatch')
    require(provenance['state_sha256'] == sha(old_state), 'original state has changed since provenance')
    require(provenance['config_sha256'] == sha(old_run / 'config.json'), 'original config has changed since provenance')
    require(provenance['checkpoint_sha256'] == sha(old_run / 'checkpoint.json'), 'original checkpoint has changed since provenance')
    require(config['resumed_from_step'] == provenance['resume_from_step'] == 7000, 'resume step mismatch')
    history = checkpoint['history']
    steps = [x['step'] for x in history]
    # The frozen resume worker appends a re-evaluation record at its restored
    # step.  It is not a second optimizer advance, but it means the terminal
    # history is not literally strictly increasing.  Reject every other
    # repeat/decrease and report this schema wart rather than silently calling
    # the history monotone.
    nonincreasing = [(i, steps[i], steps[i + 1]) for i in range(len(steps) - 1)
                     if steps[i] >= steps[i + 1]]
    require(nonincreasing == [(70, 7000, 7000)], 'unexpected non-increasing checkpoint history')
    require(history[71].get('event') == 'resume_initial_evaluation', 'controlled duplicate lacks resume event')
    require(steps[-1] == config['steps'] == 12000, 'completion does not reach requested final step')
    state = torch.load(run / 'state.pt', map_location='cpu', weights_only=True)
    require(state['step'] == 12000 and state['prefix_epoch'] == 5, 'terminal torch state mismatch')
    require(state['optimizer']['state'][0]['step'].item() == 2001, 'terminal Adam counter mismatch')
    expected = [run / f'restart_{i:07d}.json' for i in range(2000, 12000, 2000)]
    require(sorted(run.glob('restart_*.json')) == expected, 'terminal restart chain mismatch')
    require(release['terminal_status'] == 'complete' and not release['timed_out'], 'release is not a normal completion')
    require(release['cleanup_called'] and release['process_will_exit'], 'release cleanup fields missing')
    out = {
        'status': 'passed_d14_resume_input_binding_and_terminal_audit_with_documented_history_duplicate',
        'claim_boundary': 'Binds one finite resumed floating-point job and its exact saved-row Fraction audit; it is neither an exact candidate certificate nor a global exclusion.',
        'resume_step': 7000,
        'terminal_step': 12000,
        'fraction_audit_sha256': sha(args.fraction_audit),
        'resume_worker_sha256': sha(args.resume_worker),
        'original_worker_sha256': sha(args.original_worker),
        'old_state_sha256': sha(old_state),
        'terminal_state_sha256': sha(run / 'state.pt'),
        'restart_files': [x.name for x in expected],
        'stored_fraction_rows': fraction['counts']['total_saved_rows'],
        'terminal_residual_float64': history[-1]['best_max_residual_float64'],
        'history_monotonicity': {'strictly_increasing': False,
                                'sole_controlled_duplicate': {'previous_index': 70, 'resume_index': 71, 'step': 7000}},
        'release': {'cleanup_called': True, 'process_will_exit': True, 'timed_out': False},
    }
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(json.dumps(out, sort_keys=True))


if __name__ == '__main__':
    main()
