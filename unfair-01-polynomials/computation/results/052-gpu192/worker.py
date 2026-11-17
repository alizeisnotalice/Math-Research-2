#!/usr/bin/env python3
"""MPS numerical search with exhaustive legal first-fractional prefixes.

The prefix list is exact standard-library data. The MPS/CPU residuals are
numerical exploration only and never constitute an exclusion or certificate.
"""
from __future__ import annotations
import argparse, gc, hashlib, json, math, os, platform, sys, time
from pathlib import Path
import numpy as np
import torch
PREFIX_DIR = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/d14_round037/worker')
sys.path.insert(0, str(PREFIX_DIR))
import fractional_prefixes as prefixes

FROZEN_SOURCE = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/d14_round037/worker/gpu_factor_search_prefix.py')
FROZEN_SOURCE_SHA256 = 'f0044bf5a0667b1825d1b58662243bc6c9b49682c809d045680ef5c53d348d9f'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def atom_json(path: Path, value: object) -> None:
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def atom_torch(path: Path, value: object) -> None:
    tmp = path.with_suffix(path.suffix + '.tmp')
    torch.save(value, tmp)
    tmp.replace(path)


def json_parameters(a: argparse.Namespace) -> dict[str, object]:
    """Stable, JSON-safe command-line record (in particular Path -> str)."""
    return {key: (str(value) if isinstance(value, Path) else value)
            for key, value in vars(a).items()}


def load_bound_plan(path: Path, d: int, m: int, n: int, margin: float) -> tuple[list[str], dict[str, object]]:
    """Consume the independently replayed d15 survivor list in saved order."""
    plan = json.loads(path.read_text())
    assert plan.get('schema') == 'd15-necessary-prefix-plan-v1'
    p = plan.get('parameters', {})
    assert p.get('d') == d == 15 and 1 <= d < min(m, n) and p.get('margin') == margin
    codes = plan.get('survivor_codes')
    assert isinstance(codes, list) and len(codes) == plan.get('survivor_count')
    assert digest(codes) == plan.get('survivor_codes_sha256')
    assert p.get('coverage_plan_starts', 0) >= len(codes)
    for code in codes:
        assert isinstance(code, str) and len(code) == d-1 and set(code) <= {'0', '1', '2'}
        assert prefixes.legal_code(code, d)
        assert not prefixes.strict_extension_filter_reasons(code, d)
    return codes, plan


def prefix_rows(codes: list[str], starts: int, offset: int) -> tuple[np.ndarray, np.ndarray, list[str]]:
    assigned = [codes[(row + offset) % len(codes)] for row in range(starts)]
    aa, bb = zip(*(prefixes.decode(code) for code in assigned))
    return np.asarray([row[1:] for row in aa], dtype=np.float32), np.asarray([row[1:] for row in bb], dtype=np.float32), assigned


def cpu_prefix_audit(d: int, codes: list[str], plan: dict[str, object], plan_path: Path) -> dict[str, object]:
    for code in codes:
        assert prefixes.legal_code(code, d)
        assert not prefixes.strict_extension_filter_reasons(code, d)
        a, b = prefixes.decode(code)
        low = prefixes.convolution_prefix(a, b, d-1)
        assert set(low) <= {0, 1}
        assert sum(a[i]*b[d-i] for i in range(1, d)) == 0
    return {'d': d, 'raw_ternary_code_count': plan['raw_ternary_count'],
            'valid_low_prefix_count': plan['valid_low_prefix_count'],
            'valid_low_prefix_codes_sha256': plan['valid_low_prefix_codes_sha256'],
            'strict_filter_rejected_by_reason_counts': plan['strict_filter_rejected_by_reason_counts'],
            'survivor_count': len(codes), 'survivor_codes_sha256': digest(codes),
            'all_consumed_survivors_cpu_audited': True, 'float_arithmetic_used': False,
            'plan_path': str(plan_path.resolve()), 'plan_sha256': sha(plan_path),
            'plan_claim_boundary': plan['claim_boundary']}


def case_prefix_audit(av: np.ndarray, bv: np.ndarray, d: int, assigned: str) -> dict[str, object]:
    expected_a, expected_b = prefixes.decode(assigned)
    derived = []
    exact_pairs = {(0, 0): '0', (1, 0): '1', (0, 1): '2'}
    for i in range(1, d):
        pair = (av[i].item(), bv[i].item())
        assert pair in exact_pairs
        derived.append(exact_pairs[pair])
    assert av[0] == bv[0] == 1
    assert np.array_equal(av[1:d], np.asarray(expected_a[1:], dtype=av.dtype))
    assert np.array_equal(bv[1:d], np.asarray(expected_b[1:], dtype=bv.dtype))
    low = np.convolve(av[:d], bv[:d])[:d]
    cross = float(sum(av[i]*bv[d-i] for i in range(1, d)))
    assert all(value in (0.0, 1.0) for value in low) and cross == 0.0
    assert ''.join(derived) == assigned
    return {'assigned_prefix_code': assigned, 'derived_prefix_code': ''.join(derived),
            'prefix_low_convolution': low.tolist(), 'prefix_cross_at_d': cross,
            'prefix_validated_from_coefficients': True}


def cpu_tree(value):
    """Make a portable CPU-only checkpoint without trusting device storage."""
    if isinstance(value, torch.Tensor):
        return value.detach().cpu().clone()
    if isinstance(value, dict):
        return {key: cpu_tree(item) for key, item in value.items()}
    if isinstance(value, list):
        return [cpu_tree(item) for item in value]
    if isinstance(value, tuple):
        return tuple(cpu_tree(item) for item in value)
    return value


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=int, default=18); ap.add_argument('--n', type=int, default=30)
    ap.add_argument('--starts', type=int, default=16384); ap.add_argument('--steps', type=int, default=2000)
    ap.add_argument('--fractional-index', type=int, default=1); ap.add_argument('--margin', type=float, default=.1)
    ap.add_argument('--restart-every', type=int, default=0); ap.add_argument('--seed', type=int, default=20261115)
    ap.add_argument('--out', type=Path, required=True); ap.add_argument('--cpu-audit-only', action='store_true')
    ap.add_argument('--prefix-plan', type=Path, required=True)
    ap.add_argument('--seconds', type=float, default=14400)
    ap.add_argument('--allow-incomplete-coverage-smoke', action='store_true')
    a = ap.parse_args()
    assert Path('/Volumes/data').is_mount(); a.out.mkdir(parents=True, exist_ok=True)
    d, m, n, B = a.fractional_index, a.m, a.n, a.starts
    assert a.steps > 0 and 4096 <= B <= 8192 and d == 15 and 1 <= d < min(m, n) and 0 < a.margin < .5 and a.seconds > 0
    codes, plan = load_bound_plan(a.prefix_plan, d, m, n, a.margin)
    coverage_epochs = math.ceil(len(codes) / B)
    coverage_offsets = [(epoch * B) % len(codes) for epoch in range(coverage_epochs)]
    assert a.restart_every > 0
    if not a.allow_incomplete_coverage_smoke:
        assert a.steps >= coverage_epochs * a.restart_every, 'steps/restart schedule does not cover every planned window'
    audit = cpu_prefix_audit(d, codes, plan, a.prefix_plan)
    audit.update({'starts': B, 'coverage_by_start_rows': B >= len(codes),
                  'coverage_mode': 'disjoint_cyclic_windows_by_restart_epoch',
                  'coverage_epochs_required': coverage_epochs, 'coverage_epoch_offsets': coverage_offsets,
                  'full_plan_covered_after_epochs': coverage_epochs,
                  'coverage_schedule_complete': a.steps >= coverage_epochs * a.restart_every,
                  'degree_lift': {'plan_m': plan['parameters']['m'], 'plan_n': plan['parameters']['n'], 'run_m': m, 'run_n': n,
                                  'justification': 'Low-prefix and d-cross checks only use indices at most d.'},
                  'parameters': json_parameters(a),
                  'worker_sha256': sha(Path(__file__)), 'prefix_module_sha256': sha(Path(prefixes.__file__)),
                  'frozen_source_path': str(FROZEN_SOURCE), 'frozen_source_sha256': FROZEN_SOURCE_SHA256})
    atom_json(a.out/'cpu-prefix-audit.json', audit)
    if a.cpu_audit_only:
        print(json.dumps(audit, indent=2)); return
    assert torch.backends.mps.is_available()
    torch.set_num_threads(1); torch.manual_seed(a.seed); dev = torch.device('mps')
    free = torch.nn.Parameter(torch.rand((B, m+n-2*d-2), device=dev)*.3)
    theta = torch.nn.Parameter(torch.randn((B, 1), device=dev))
    one = torch.ones((B, 1), device=dev)
    incidence_cpu = np.zeros(((m+1)*(n+1), m+n+1), dtype=np.float32)
    for i in range(m+1):
        for j in range(n+1): incidence_cpu[i*(n+1)+j, i+j] = 1
    incidence = torch.tensor(incidence_cpu, device=dev)
    opt = torch.optim.Adam([free, theta], lr=.035)
    metadata = {**json_parameters(a), 'torch': torch.__version__, 'numpy': np.__version__,
                'device': 'mps', 'platform': platform.platform(), 'dtype': 'float32', 'pid': os.getpid(),
                'script_sha256': sha(Path(__file__)), 'prefix_module_sha256': sha(Path(prefixes.__file__)),
                'prefix_d': d, 'prefix_survivor_count': len(codes), 'prefix_codes_sha256': audit['survivor_codes_sha256'],
                'frozen_source_path': str(FROZEN_SOURCE), 'frozen_source_sha256': FROZEN_SOURCE_SHA256,
                'prefix_plan_path': str(a.prefix_plan.resolve()), 'prefix_plan_sha256': sha(a.prefix_plan),
                'coverage_mode': 'disjoint_cyclic_windows_by_restart_epoch',
                'coverage_epochs_required': coverage_epochs, 'coverage_epoch_offsets': coverage_offsets,
                'coverage_schedule_complete': a.steps >= coverage_epochs * a.restart_every,
                'degree_lift': {'plan_m': plan['parameters']['m'], 'plan_n': plan['parameters']['n'], 'run_m': m, 'run_n': n,
                                'justification': 'Low-prefix and d-cross checks only use indices at most d.'},
                'scope': 'Filtered legal low prefixes; projected numerical MPS search only.'}
    atom_json(a.out/'config.json', metadata)
    assigned: list[str] = []; pa = pb = None; cases: list[dict[str, object]] = []; history = []; window_candidates = []
    start = time.monotonic(); restart_epoch = 0

    def rotate_prefixes(epoch: int) -> None:
        nonlocal assigned, pa, pb
        pa_cpu, pb_cpu, assigned = prefix_rows(codes, B, epoch * B)
        pa = torch.tensor(pa_cpu, device=dev); pb = torch.tensor(pb_cpu, device=dev)

    def factors():
        assert pa is not None and pb is not None
        t = a.margin + (1-2*a.margin)*torch.sigmoid(theta)
        aa = torch.cat((one, pa, t, free[:, :m-d-1], one), dim=1)
        bb = torch.cat((one, pb, 1-t, free[:, m-d-1:], one), dim=1)
        return aa, bb

    def save_state(step: int) -> None:
        # This is a local, resumable temporary artifact, not a mathematical certificate.
        atom_torch(a.out/'state.pt', {'step': step, 'free': free.detach().cpu().clone(),
                    'theta': theta.detach().cpu().clone(), 'optimizer': cpu_tree(opt.state_dict()),
                    'prefix_epoch': restart_epoch, 'prefix_offset_mod_count': (restart_epoch * B) % len(codes),
                    'prefix_codes_sha256': audit['survivor_codes_sha256'],
                    'parameters': json_parameters(a)})

    rotate_prefixes(restart_epoch)
    last_completed_step = -1; timed_out = False; terminal_status = 'running'
    try:
        for step in range(a.steps+1):
            if time.monotonic() - start >= a.seconds:
                timed_out = True; terminal_status = 'partial_time_limit'; break
            if a.restart_every and 0 < step < a.steps and step % a.restart_every == 0:
                # Preserve the old cohort *before* its buffers/optimizer are reset.
                window_record = {'window_epoch': restart_epoch,
                           'window_offset_mod_count': (restart_epoch * B) % len(codes),
                           'closing_step': history[-1]['step'] if history else None,
                           'candidates': cases}
                window_candidates.append(window_record)
                atom_json(a.out/'window-candidates.json', {'schema': 'd15-window-candidates-v1',
                          'seed': a.seed, 'starts': B, 'restart_every': a.restart_every,
                          'prefix_codes_sha256': audit['survivor_codes_sha256'], 'windows': window_candidates,
                          'claim_boundary': 'Saved floating-point candidate rows only; no exact certificate.'})
                atom_json(a.out/f'restart_{step:07d}.json', {
                    'old_prefix_epoch': restart_epoch,
                    'old_prefix_offset_mod_count': (restart_epoch * B) % len(codes),
                    'next_prefix_epoch': restart_epoch + 1,
                    'next_prefix_offset_mod_count': ((restart_epoch + 1) * B) % len(codes),
                    'history': history[-1:], 'cases': cases})
                restart_epoch += 1; rotate_prefixes(restart_epoch)
                with torch.no_grad(): free.uniform_(0, .3); theta.normal_()
                opt.state.clear()
            opt.zero_grad(set_to_none=True); aa, bb = factors()
            conv = (aa[:, :, None]*bb[:, None, :]).reshape(B, -1) @ incidence
            target = (conv.detach() >= .5).to(conv.dtype); residual = conv-target; per = (residual*residual).mean(dim=1)
            loss = per.sum(); loss.backward(); opt.step()
            with torch.no_grad(): free.clamp_(0, 1)
            last_completed_step = step
            if step % 100 == 0 or step == a.steps:
                torch.mps.synchronize()
                with torch.no_grad():
                    scores, ix = torch.topk(per, 16, largest=False); aa, bb = factors()
                    va, vb = aa[ix].detach().cpu().numpy(), bb[ix].detach().cpu().numpy(); cases = []
                    for row, (av, bv) in enumerate(zip(va, vb)):
                        source_index = int(ix[row].cpu()); prefix_case = case_prefix_audit(av, bv, d, assigned[source_index])
                        cv = np.convolve(av.astype(np.float64), bv.astype(np.float64)); tar = (cv >= .5).astype(np.float64)
                        cases.append({**prefix_case, 'A': av.tolist(), 'B': bv.tolist(),
                                      'max_residual_float64': float(np.max(np.abs(cv-tar))),
                                      'l2_residual_float64': float(np.linalg.norm(cv-tar)), 'status': 'numerical_only'})
                    record = {'step': step, 'elapsed_seconds': time.monotonic()-start, 'restart_epoch': restart_epoch,
                              'min_mse_gpu': float(scores[0].cpu()),
                              'best_max_residual_float64': min(c['max_residual_float64'] for c in cases),
                              'allocated_bytes': torch.mps.current_allocated_memory(), 'driver_bytes': torch.mps.driver_allocated_memory()}
                    history.append(record); atom_json(a.out/'checkpoint.json', {'history': history, 'cases': cases})
                    if step % 500 == 0 or step == a.steps:
                        save_state(step)
                    print(json.dumps(record), flush=True)
        if timed_out:
            save_state(last_completed_step)
            terminal_status = 'partial_time_limit'
        else:
            terminal_status = 'complete'
        if cases:
            window_candidates.append({'window_epoch': restart_epoch,
                                      'window_offset_mod_count': (restart_epoch * B) % len(codes),
                                      'closing_step': history[-1]['step'] if history else None,
                                      'candidates': cases})
        atom_json(a.out/'window-candidates.json', {'schema': 'd15-window-candidates-v1',
                  'seed': a.seed, 'starts': B, 'restart_every': a.restart_every,
                  'prefix_codes_sha256': audit['survivor_codes_sha256'], 'windows': window_candidates,
                  'claim_boundary': 'Saved floating-point candidate rows only; no exact certificate.'})
        atom_json(a.out/'result.json', {'complete': not timed_out, 'status': terminal_status,
            'last_completed_step': last_completed_step, 'history': history, 'cases': cases,
            'note': 'Floating-point search is neither an exclusion certificate nor an exact counterexample.'})
    except BaseException as exc:
        terminal_status = 'failed_exception'
        atom_json(a.out/'result.json', {'complete': False, 'status': terminal_status,
            'last_completed_step': last_completed_step, 'history': history, 'cases': cases,
            'exception_type': type(exc).__name__, 'note': 'Numerical batch failed and is not complete.'})
        raise
    finally:
        del opt, free, theta, one, incidence, pa, pb
        gc.collect(); torch.mps.synchronize(); torch.mps.empty_cache()
        atom_json(a.out/'release.json', {'cleanup_called': True, 'process_will_exit': True,
            'terminal_status': terminal_status, 'timed_out': timed_out, 'last_completed_step': last_completed_step,
            'remaining_tensor_bytes_before_exit': torch.mps.current_allocated_memory(), 'pid': os.getpid()})


if __name__ == '__main__': main()
