#!/usr/bin/env python3
"""Resume an interrupted audited-prefix MPS exploration from a CPU checkpoint.

This is deliberately a new worker.  It never mutates the interrupted source
directory: it validates its saved state, copies only pre-state restart records,
and records which later observations were discarded.  The computation remains
floating-point exploration, not an exclusion or an exact certificate.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import math
import os
import platform
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
SOURCE = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/d15_prefix_plan_39/gpu_factor_search_d15.py')
SOURCE_SHA = '5485fd6233a3dd62f1936dd20c23d2ab4effd9bec5187e5958cbcd8c91395de0'
PLAN_SHA = 'aa1d442efdb8efe1c1bd65a99f7ab8021222553008ecc915b13f9800323e5a7a'
PLAN_CODES_SHA = '80f10e03d019aad3513165f940e45c311c50c3fb7dd981c3a4ab62cad0e009df'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def durable_replace(tmp: Path, path: Path) -> None:
    fd = os.open(tmp, os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)
    os.replace(tmp, path)
    directory = os.open(path.parent, os.O_RDONLY)
    try: os.fsync(directory)
    finally: os.close(directory)


def atom(path: Path, obj: object) -> None:
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')
    durable_replace(tmp, path)


def atom_torch(path: Path, obj: object) -> None:
    """Write a checkpoint atomically, so a second reboot cannot half-write it."""
    tmp = path.with_suffix(path.suffix + '.tmp')
    torch.save(obj, tmp)
    durable_replace(tmp, path)


def load_source():
    if not SOURCE.is_file() or sha(SOURCE) != SOURCE_SHA:
        raise ValueError('recorded d15 source is absent or has a different SHA-256')
    spec = importlib.util.spec_from_file_location('frozen_d15', SOURCE)
    if spec is None or spec.loader is None:
        raise ValueError('could not import recorded d14 worker')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def params(a: argparse.Namespace) -> dict[str, object]:
    return {'m': a.m, 'n': a.n, 'starts': a.starts, 'steps': a.steps,
            'fractional_index': a.fractional_index, 'margin': a.margin,
            'restart_every': a.restart_every, 'seed': a.seed}


def require(condition: bool, text: str) -> None:
    if not condition:
        raise ValueError(text)


def finite_cpu_tensor(x: object, shape: tuple[int, ...], label: str) -> None:
    require(isinstance(x, torch.Tensor), f'{label} is not a tensor')
    require(x.device.type == 'cpu', f'{label} is not stored on CPU')
    require(tuple(x.shape) == shape, f'{label} has shape {tuple(x.shape)}, expected {shape}')
    require(x.dtype == torch.float32, f'{label} has dtype {x.dtype}, expected float32')
    require(bool(torch.isfinite(x).all()), f'{label} has non-finite entries')


def scalar_step(value: object, label: str) -> int:
    if isinstance(value, torch.Tensor):
        require(value.device.type == 'cpu' and value.numel() == 1 and bool(torch.isfinite(value).all()), f'{label} invalid')
        value = value.item()
    require(isinstance(value, (int, float)) and math.isfinite(value) and int(value) == value and value >= 0, f'{label} invalid')
    return int(value)


def validate_optimizer(state: dict, free_shape: tuple[int, ...], theta_shape: tuple[int, ...], expected_step: int) -> None:
    require(isinstance(state, dict) and set(state) == {'state', 'param_groups'}, 'optimizer state schema differs')
    require(isinstance(state['param_groups'], list) and len(state['param_groups']) == 1, 'optimizer has not exactly one parameter group')
    group = state['param_groups'][0]
    require(group.get('params') == [0, 1], 'optimizer parameter identifiers differ')
    require(group.get('lr') == .035 and group.get('betas') == (0.9, 0.999) and group.get('eps') == 1e-8, 'Adam learning parameters differ')
    require(group.get('weight_decay') == 0 and group.get('amsgrad') is False and group.get('maximize') is False, 'Adam regularization/mode differs')
    require(group.get('foreach') is None and group.get('capturable') is False and group.get('differentiable') is False and group.get('fused') is None and group.get('decoupled_weight_decay') is False, 'Adam implementation options differ')
    slots = state['state']
    require(set(slots) == {0, 1}, 'optimizer slots do not match free/theta')
    for key, shape in [(0, free_shape), (1, theta_shape)]:
        slot = slots[key]
        require(set(slot) == {'step', 'exp_avg', 'exp_avg_sq'}, f'optimizer slot {key} schema differs')
        finite_cpu_tensor(slot['exp_avg'], shape, f'optimizer[{key}].exp_avg')
        finite_cpu_tensor(slot['exp_avg_sq'], shape, f'optimizer[{key}].exp_avg_sq')
        require(scalar_step(slot['step'], f'optimizer[{key}].step') == expected_step, f'optimizer[{key}].step does not match saved state step')


def load_and_validate_resume(a: argparse.Namespace, old, old_state: Path) -> dict[str, object]:
    run = a.resume_run.resolve()
    require(run.is_dir(), 'resume-run is not a directory')
    require(old_state.is_file(), 'resume-state is not a file')
    config_path = run / 'config.json'
    checkpoint_path = run / 'checkpoint.json'
    require(config_path.is_file() and checkpoint_path.is_file(), 'resume-run lacks config.json or checkpoint.json')
    old_config = json.loads(config_path.read_text())
    current = params(a)
    for key, value in current.items():
        require(old_config.get(key) == value, f'resume configuration mismatch for {key}: {old_config.get(key)!r} != {value!r}')
    require(old_config.get('script_sha256') == SOURCE_SHA, 'resume config does not bind the recorded d15 source')
    require(old_config.get('prefix_plan_sha256') == PLAN_SHA, 'resume config does not bind the recorded d15 plan')
    require(Path(old_config.get('prefix_plan_path', '')).resolve() == a.prefix_plan.resolve(), 'resume config plan path mismatch')
    require(old_config.get('prefix_module_sha256') == sha(Path(old.prefixes.__file__)), 'resume prefix module hash mismatch')
    data = torch.load(old_state, map_location='cpu', weights_only=True)
    required = {'step', 'free', 'theta', 'optimizer', 'prefix_epoch', 'prefix_offset_mod_count', 'prefix_codes_sha256', 'parameters'}
    require(set(data) == required, f'resume state keys differ: {set(data)}')
    for key, value in current.items():
        require(data['parameters'].get(key) == value, f'state parameters mismatch for {key}')
    codes, plan = old.load_bound_plan(a.prefix_plan, a.fractional_index, a.m, a.n, a.margin)
    digest = hashlib.sha256(json.dumps(codes, separators=(',', ':')).encode()).hexdigest()
    require(digest == PLAN_CODES_SHA, 'bound d15 plan code order differs from the audited list')
    require(data['prefix_codes_sha256'] == digest == old_config.get('prefix_codes_sha256'), 'prefix code SHA mismatch')
    step = data['step']
    require(isinstance(step, int) and 0 <= step < a.steps, 'state step must be an integer in [0, steps)')
    epoch = data['prefix_epoch']
    require(isinstance(epoch, int) and epoch >= 0 and data['prefix_offset_mod_count'] == epoch % len(codes), 'prefix epoch/offset is inconsistent')
    expected_epoch = step // a.restart_every if a.restart_every else 0
    require(epoch == expected_epoch, f'checkpoint epoch {epoch} inconsistent with step {step}')
    free_shape = (a.starts, a.m + a.n - 2 * a.fractional_index - 2)
    theta_shape = (a.starts, 1)
    finite_cpu_tensor(data['free'], free_shape, 'free')
    require(bool(((data['free'] >= 0) & (data['free'] <= 1)).all()), 'free is outside the original clamp interval [0,1]')
    finite_cpu_tensor(data['theta'], theta_shape, 'theta')
    expected_adam_step = (step % a.restart_every) + 1 if a.restart_every else step + 1
    validate_optimizer(data['optimizer'], free_shape, theta_shape, expected_adam_step)
    checkpoint = json.loads(checkpoint_path.read_text())
    history = checkpoint.get('history')
    require(isinstance(history, list), 'checkpoint history is absent')
    require(all(isinstance(x, dict) and isinstance(x.get('step'), int) for x in history), 'checkpoint history has a malformed step')
    history_steps = [x['step'] for x in history]
    require(all(left < right for left, right in zip(history_steps, history_steps[1:])), 'checkpoint history steps are not strictly increasing')
    kept = [x for x in history if isinstance(x, dict) and x.get('step', -1) <= step]
    discarded = [x for x in history if not (isinstance(x, dict) and x.get('step', -1) <= step)]
    require(kept and kept[-1]['step'] == step, 'checkpoint has no history record at the saved state step')
    require(all(x.get('step', -1) <= step for x in kept), 'history truncation failed')
    expected_restart_steps = list(range(a.restart_every, step + 1, a.restart_every)) if a.restart_every else []
    actual_restart_paths = sorted(run.glob('restart_*.json'))
    actual_restart_steps = []
    for path in actual_restart_paths:
        try:
            restart_step = int(path.stem.split('_')[1])
        except (IndexError, ValueError) as exc:
            raise ValueError(f'malformed restart filename {path.name}') from exc
        actual_restart_steps.append(restart_step)
        if restart_step > step:
            continue
        doc = json.loads(path.read_text())
        require(set(doc) == {'old_prefix_epoch', 'old_prefix_offset_mod_count', 'next_prefix_epoch', 'next_prefix_offset_mod_count', 'history', 'cases'}, f'{path.name} schema differs')
        old_epoch = restart_step // a.restart_every - 1
        require(doc['old_prefix_epoch'] == old_epoch and doc['next_prefix_epoch'] == old_epoch + 1, f'{path.name} epoch chain differs')
        require(doc['old_prefix_offset_mod_count'] == old_epoch % len(codes) and doc['next_prefix_offset_mod_count'] == (old_epoch + 1) % len(codes), f'{path.name} prefix offset chain differs')
        require(isinstance(doc['history'], list) and isinstance(doc['cases'], list), f'{path.name} history/cases malformed')
    require([x for x in actual_restart_steps if x <= step] == expected_restart_steps, 'restart record chain is incomplete or has extra records')
    restart_paths = [run / f'restart_{restart_step:07d}.json' for restart_step in expected_restart_steps]
    # Rebuild all factors on CPU for the restored epoch.  This is an exact
    # prefix-order check plus endpoint/complement check, not a GPU evaluation.
    pa, pb, assigned = old.prefix_rows(codes, a.starts, epoch)
    t = a.margin + (1 - 2 * a.margin) / (1 + np.exp(-data['theta'].numpy()))
    factor_a = np.concatenate((np.ones((a.starts, 1), dtype=np.float32), pa, t.astype(np.float32),
                               data['free'].numpy()[:, :a.m-a.fractional_index-1], np.ones((a.starts, 1), dtype=np.float32)), axis=1)
    factor_b = np.concatenate((np.ones((a.starts, 1), dtype=np.float32), pb, (1-t).astype(np.float32),
                               data['free'].numpy()[:, a.m-a.fractional_index-1:], np.ones((a.starts, 1), dtype=np.float32)), axis=1)
    require(factor_a.shape == (a.starts, a.m + 1) and factor_b.shape == (a.starts, a.n + 1), 'rebuilt factor shapes differ')
    require(np.array_equal(factor_a[:, 1:a.fractional_index], pa) and np.array_equal(factor_b[:, 1:a.fractional_index], pb), 'rebuilt factors lost saved prefix order')
    require(np.all(factor_a[:, 0] == 1) and np.all(factor_a[:, -1] == 1) and np.all(factor_b[:, 0] == 1) and np.all(factor_b[:, -1] == 1), 'rebuilt factor endpoints differ')
    require(np.allclose(factor_a[:, a.fractional_index] + factor_b[:, a.fractional_index], 1, rtol=0, atol=2e-7), 'fractional complement differs')
    return {'state': data, 'codes': codes, 'plan': plan, 'assigned': assigned, 'step': step, 'epoch': epoch, 'history': kept,
            'discarded_history': discarded, 'restart_paths': restart_paths,
            'old_config': old_config, 'state_sha256': sha(old_state),
            'checkpoint_sha256': sha(checkpoint_path), 'config_sha256': sha(config_path)}


def write_resume_provenance(a: argparse.Namespace, record: dict[str, object], out: Path) -> None:
    copied = []
    for source in record['restart_paths']:
        target = out / source.name
        shutil.copy2(source, target)
        copied.append({'name': source.name, 'sha256': sha(target)})
    atom(out / 'resume-provenance.json', {
        'status': 'validated_resume_state', 'resume_run': str(a.resume_run.resolve()),
        'resume_state': str(a.resume_state.resolve()), 'source_worker': str(SOURCE),
        'source_worker_sha256': SOURCE_SHA, 'source_plan': str(a.prefix_plan.resolve()),
        'source_plan_sha256': PLAN_SHA, 'source_plan_codes_sha256': PLAN_CODES_SHA,
        'resume_worker_sha256': sha(Path(__file__)),
        'state_sha256': record['state_sha256'], 'config_sha256': record['config_sha256'],
        'checkpoint_sha256': record['checkpoint_sha256'], 'resume_from_step': record['step'],
        'next_step': record['step'] + 1, 'prefix_epoch': record['epoch'],
        'kept_history_steps': [x['step'] for x in record['history']],
        'discarded_history_steps': [x.get('step') for x in record['discarded_history']],
        'copied_restart_records': copied,
        'recomputed_initial_cases': True,
        'rng_note': 'The source worker did not checkpoint RNG state. Adam parameters/moments are restored, but a later restart uses a new deterministic seed sequence and is not bitwise equivalent to the interrupted process; no RNG state is invented.',
        'claim_boundary': 'Numerical resume provenance only; no factorization, exclusion, or exact counterexample is certified.'
    })


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=int, required=True); ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--starts', type=int, required=True); ap.add_argument('--steps', type=int, required=True)
    ap.add_argument('--fractional-index', type=int, required=True); ap.add_argument('--margin', type=float, required=True)
    ap.add_argument('--restart-every', type=int, required=True); ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--out', type=Path, required=True); ap.add_argument('--resume-run', type=Path, required=True)
    ap.add_argument('--resume-state', type=Path); ap.add_argument('--validate-resume-only', action='store_true')
    ap.add_argument('--prefix-plan', type=Path, required=True)
    ap.add_argument('--preflight-record', type=Path, help='write the CPU-only preflight receipt before exiting')
    ap.add_argument('--seconds', type=float, default=14400, help='hard wall-clock bound after restore')
    a = ap.parse_args()
    require(Path('/Volumes/data').is_mount(), 'data is not mounted')
    require(a.steps > 0 and a.starts >= 16 and a.fractional_index == 15 and 15 < min(a.m, a.n) and 0 < a.margin < .5 and a.seconds > 0, 'invalid search parameters')
    if a.resume_state is None:
        a.resume_state = a.resume_run / 'state.pt'
    old = load_source()
    record = load_and_validate_resume(a, old, a.resume_state)
    summary = {'status': 'resume_preflight_passed', 'step': record['step'], 'next_step': record['step'] + 1,
               'epoch': record['epoch'], 'kept_history': len(record['history']),
               'discarded_history_steps': [x.get('step') for x in record['discarded_history']],
               'state_sha256': record['state_sha256'], 'source_worker_sha256': SOURCE_SHA,
               'source_plan_sha256': PLAN_SHA, 'source_plan_codes_sha256': PLAN_CODES_SHA,
               'resume_worker_sha256': sha(Path(__file__)), 'torch': torch.__version__}
    if a.preflight_record is not None:
        require(a.preflight_record.parent.is_dir(), 'preflight-record parent does not exist')
        atom(a.preflight_record, summary)
    if a.validate_resume_only:
        print(json.dumps(summary, sort_keys=True)); return
    require(not a.out.exists() or not any(a.out.iterdir()), 'out exists and is nonempty')
    a.out.mkdir(parents=True, exist_ok=True)
    write_resume_provenance(a, record, a.out)
    require(torch.backends.mps.is_available(), 'MPS is unavailable')
    torch.set_num_threads(1); torch.manual_seed(a.seed); dev = torch.device('mps')
    codes = record['codes']; d = a.fractional_index
    audit = old.cpu_prefix_audit(d, codes, record['plan'], a.prefix_plan)
    audit.update({'starts': a.starts, 'coverage_by_start_rows': True, 'parameters': params(a),
                  'worker_sha256': sha(Path(__file__)), 'original_worker_sha256': SOURCE_SHA,
                  'prefix_module_sha256': sha(Path(old.prefixes.__file__)), 'resumed_from_step': record['step']})
    atom(a.out / 'cpu-prefix-audit.json', audit)
    atom(a.out / 'config.json', {**params(a), 'out': str(a.out), 'resume_run': str(a.resume_run.resolve()),
                                  'resume_state': str(a.resume_state.resolve()), 'seconds_limit': a.seconds,
                                  'torch': torch.__version__, 'numpy': np.__version__, 'device': 'mps',
                                  'platform': platform.platform(), 'dtype': 'float32', 'pid': os.getpid(),
                                  'script_sha256': sha(Path(__file__)), 'original_worker_sha256': SOURCE_SHA,
                                  'prefix_module_sha256': sha(Path(old.prefixes.__file__)), 'prefix_d': d,
                                  'prefix_survivor_count': len(codes), 'prefix_codes_sha256': audit['survivor_codes_sha256'],
                                  'prefix_plan_path': str(a.prefix_plan.resolve()), 'prefix_plan_sha256': PLAN_SHA,
                                  'resumed_from_step': record['step'], 'scope': 'Resumed filtered-prefix numerical MPS search only.'})
    free = torch.nn.Parameter(record['state']['free'].to(dev)); theta = torch.nn.Parameter(record['state']['theta'].to(dev))
    one = torch.ones((a.starts, 1), device=dev)
    incidence_cpu = np.zeros(((a.m + 1) * (a.n + 1), a.m + a.n + 1), dtype=np.float32)
    for i in range(a.m + 1):
        for j in range(a.n + 1): incidence_cpu[i * (a.n + 1) + j, i + j] = 1
    incidence = torch.tensor(incidence_cpu, device=dev)
    opt = torch.optim.Adam([free, theta], lr=.035); opt.load_state_dict(record['state']['optimizer'])
    pa = pb = None; assigned: list[str] = []; epoch = record['epoch']; history = list(record['history']); cases: list[dict[str, object]] = []

    def rotate(which: int) -> None:
        nonlocal pa, pb, assigned
        aa, bb, assigned = old.prefix_rows(codes, a.starts, which)
        pa, pb = torch.tensor(aa, device=dev), torch.tensor(bb, device=dev)

    def factors():
        t = a.margin + (1 - 2 * a.margin) * torch.sigmoid(theta)
        aa = torch.cat((one, pa, t, free[:, :a.m-d-1], one), dim=1)
        bb = torch.cat((one, pb, 1-t, free[:, a.m-d-1:], one), dim=1)
        return aa, bb

    def save_state(step: int) -> None:
        atom_torch(a.out / 'state.pt', {'step': step, 'free': free.detach().cpu().clone(),
            'theta': theta.detach().cpu().clone(), 'optimizer': old.cpu_tree(opt.state_dict()),
            'prefix_epoch': epoch, 'prefix_offset_mod_count': epoch % len(codes),
            'prefix_codes_sha256': audit['survivor_codes_sha256'], 'parameters': params(a)})

    def evaluate(step: int, tag: str) -> None:
        nonlocal cases
        with torch.no_grad():
            aa, bb = factors(); conv = (aa[:, :, None] * bb[:, None, :]).reshape(a.starts, -1) @ incidence
            target = (conv >= .5).to(conv.dtype); per = ((conv-target)**2).mean(dim=1); scores, ix = torch.topk(per, 16, largest=False)
            va, vb = aa[ix].detach().cpu().numpy(), bb[ix].detach().cpu().numpy(); cases = []
            for row, (av, bv) in enumerate(zip(va, vb)):
                prefix = old.case_prefix_audit(av, bv, d, assigned[int(ix[row].cpu())])
                cv = np.convolve(av.astype(np.float64), bv.astype(np.float64)); target_cpu = (cv >= .5).astype(np.float64)
                cases.append({**prefix, 'A': av.tolist(), 'B': bv.tolist(), 'max_residual_float64': float(np.max(np.abs(cv-target_cpu))), 'status': 'numerical_only'})
            item = {'step': step, 'event': tag, 'restart_epoch': epoch, 'min_mse_gpu': float(scores[0].cpu()),
                    'best_max_residual_float64': min(x['max_residual_float64'] for x in cases),
                    'recomputed_from_restored_state': tag == 'resume_initial_evaluation'}
            if tag == 'resume_initial_evaluation' and history and history[-1]['step'] == step:
                history[-1] = item
            else:
                history.append(item)
            atom(a.out / 'checkpoint.json', {'history': history, 'cases': cases})

    rotate(epoch); evaluate(record['step'], 'resume_initial_evaluation'); save_state(record['step'])
    started = time.monotonic(); last_completed = record['step']; timed_out = False; terminal_status = 'running'
    try:
        for step in range(record['step'] + 1, a.steps + 1):
            if time.monotonic() - started >= a.seconds:
                timed_out = True; terminal_status = 'partial_time_limit'; break
            if a.restart_every and step % a.restart_every == 0 and step < a.steps:
                atom(a.out / f'restart_{step:07d}.json', {'old_prefix_epoch': epoch, 'next_prefix_epoch': epoch + 1, 'history': history[-1:], 'cases': cases, 'resume_generated': True})
                epoch += 1; rotate(epoch)
                with torch.no_grad(): free.uniform_(0, .3); theta.normal_()
                opt.state.clear()
            opt.zero_grad(set_to_none=True); aa, bb = factors(); conv = (aa[:, :, None] * bb[:, None, :]).reshape(a.starts, -1) @ incidence
            target = (conv.detach() >= .5).to(conv.dtype); per = ((conv-target)**2).mean(dim=1); per.sum().backward(); opt.step()
            with torch.no_grad(): free.clamp_(0, 1)
            last_completed = step
            if step % 100 == 0 or step == a.steps:
                torch.mps.synchronize(); evaluate(step, 'resumed_progress')
                print(json.dumps(history[-1]), flush=True)
            if step % 500 == 0 or step == a.steps:
                save_state(step)
        if timed_out:
            save_state(last_completed)
            terminal_status = 'partial_time_limit'
        else:
            terminal_status = 'complete'
        atom(a.out / 'result.json', {'complete': not timed_out, 'status': terminal_status,
            'last_completed_step': last_completed, 'history': history, 'cases': cases,
            'resumed_from_step': record['step'], 'note': 'Floating-point exploration only.'})
    except BaseException as exc:
        terminal_status = 'failed_exception'
        atom(a.out / 'result.json', {'complete': False, 'status': terminal_status,
            'last_completed_step': last_completed, 'history': history, 'cases': cases,
            'resumed_from_step': record['step'], 'exception_type': type(exc).__name__,
            'note': 'The numerical run failed; it is not a completed batch.'})
        raise
    finally:
        del opt, free, theta, one, incidence, pa, pb
        gc.collect(); torch.mps.synchronize(); torch.mps.empty_cache()
        atom(a.out / 'release.json', {'cleanup_called': True, 'process_will_exit': True,
            'terminal_status': terminal_status, 'timed_out': timed_out,
            'last_completed_step': last_completed, 'remaining_tensor_bytes_before_exit': torch.mps.current_allocated_memory(), 'pid': os.getpid()})


if __name__ == '__main__':
    main()
