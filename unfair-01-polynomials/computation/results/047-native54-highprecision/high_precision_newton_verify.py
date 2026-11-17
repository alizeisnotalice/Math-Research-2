#!/usr/bin/env python3
"""Independent high-precision numerical diagnostic for frozen degree-54 masks.

The target C is exact integer data; the supplied group roots are only binary64
seeds.  Each seed is refined by Newton iteration against C at arbitrary mpmath
precision.  Results are numerical diagnostics, never root-allocation or
algebraic certificates.  In the exceptional case that rounded 0/1 factors
multiply *exactly* to C using integer arithmetic, that integer product is
recorded separately as a genuine fair decomposition of C.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import mpmath as mp


def polyval_asc(c, z):
    value = mp.mpc(0)
    for a in reversed(c):
        value = value*z + a
    return value


def derivative_asc(c, z):
    value = mp.mpc(0)
    for k in range(len(c)-1, 0, -1):
        value = value*z + k*c[k]
    return value


def newton(seed, c, tol, max_steps):
    z = mp.mpc(mp.mpf(str(seed[0])), mp.mpf(str(seed[1])))
    for step in range(1, max_steps + 1):
        value = polyval_asc(c, z)
        deriv = derivative_asc(c, z)
        if not mp.isfinite(value.real) or not mp.isfinite(value.imag) or abs(deriv) == 0:
            return None, {'status': 'nonfinite_or_critical', 'steps': step}
        delta = value / deriv
        z -= delta
        if abs(delta) <= tol * max(1, abs(z)) and abs(polyval_asc(c, z)) <= tol:
            return z, {'status': 'converged', 'steps': step, 'residual': abs(polyval_asc(c, z)), 'update': abs(delta)}
    return None, {'status': 'iteration_limit', 'steps': max_steps, 'residual': abs(polyval_asc(c, z))}


def convolve(p, q):
    out = [mp.mpc(0)] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            out[i+j] += pi*qj
    return out


def group_factor(group):
    out = [mp.mpc(1)]
    for z in group:
        out = convolve(out, [mp.mpc(1), -1/z])
    return out


def bits_for_mask(groups, mask):
    if not isinstance(mask, int) or not (0 <= mask < 1 << (groups-1)):
        raise ValueError('mask outside frozen half-space')
    # Reimplemented from the documented frozen convention: group 0 is B,
    # while the remaining bits run high group index to low bit significance.
    return [0] + [(mask >> (groups-1-g)) & 1 for g in range(1, groups)]


def max_abs(values):
    return max((abs(x) for x in values), default=mp.mpf(0))


def min_max_real(values):
    real = [x.real for x in values]
    return min(real), max(real)


def distance_binary(values):
    return max((min(abs(z.real), abs(z.real-1)) for z in values), default=mp.mpf(0))


def integer_convolution(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i+j] += ai*bj
    return out


def scientific(x, digits=18):
    return mp.nstr(x, digits, min_fixed=0, max_fixed=0)


def diagnose_task(task, masks, dps, tol_exp, max_steps):
    c = [int(x) for x in task['C']]
    tol = mp.mpf(10) ** (-tol_exp)
    roots, logs = [], []
    for group in task['roots']:
        refined_group = []
        for seed in group:
            z, log = newton(seed, c, tol, max_steps)
            logs.append(log)
            if z is None:
                return {'task': task['id'], 'root_status': 'failed', 'root_logs': logs}
            refined_group.append(z)
        roots.append(refined_group)
    flat = [z for group in roots for z in group]
    separation = min(abs(flat[i]-flat[j]) for i in range(len(flat)) for j in range(i))
    conjugate_defect = mp.mpf(0)
    for group in roots:
        if len(group) == 2:
            conjugate_defect = max(conjugate_defect, abs(group[0]-mp.conj(group[1])))
        elif len(group) == 1:
            conjugate_defect = max(conjugate_defect, abs(group[0].imag))
        else:
            return {'task': task['id'], 'root_status': 'invalid_group_cardinality', 'root_logs': logs}
    factors = [group_factor(group) for group in roots]
    all_factor = [mp.mpc(1)]
    for factor in factors:
        all_factor = convolve(all_factor, factor)
    control_error = max_abs([all_factor[k]-c[k] for k in range(len(c))])
    converged_logs = [row for row in logs if row['status'] == 'converged']
    return {
        'task': task['id'], 'root_status': 'converged', 'roots': roots,
        'factors': factors, 'root_logs': logs, 'min_root_separation': separation,
        'max_conjugate_defect': conjugate_defect, 'all_groups_target_error': control_error,
        'max_newton_steps_used': max(row['steps'] for row in converged_logs),
        'max_final_root_residual': max(row['residual'] for row in converged_logs),
        'tolerance': tol,
    }


def diagnose_mask(root_data, task, mask, tol):
    if root_data['root_status'] != 'converged':
        return {'task': task['id'], 'mask': mask, 'classification': 'root_refinement_failed'}
    bits = bits_for_mask(task['groups'], mask)
    a, b = [mp.mpc(1)], [mp.mpc(1)]
    for bit, factor in zip(bits, root_data['factors']):
        if bit:
            a = convolve(a, factor)
        else:
            b = convolve(b, factor)
    ab = convolve(a, b)
    c = [int(x) for x in task['C']]
    product_error = max_abs([ab[k]-c[k] for k in range(len(c))])
    values = a+b
    imag = max_abs([z.imag for z in values] + [z.imag for z in ab])
    lo, hi = min_max_real(values)
    dbin = distance_binary(values)
    rounded_a = [int(mp.nint(z.real)) for z in a]
    rounded_b = [int(mp.nint(z.real)) for z in b]
    rounding_error = max_abs([z.real-r for z, r in zip(a, rounded_a)] + [z.real-r for z, r in zip(b, rounded_b)])
    exact_round_product = (
        all(v in (0, 1) for v in rounded_a+rounded_b)
        and integer_convolution(rounded_a, rounded_b) == c
    )
    if imag > tol or product_error > tol:
        label = 'refinement_inconsistent'
    elif lo < -tol:
        label = 'negative'
    elif hi > 1+tol:
        label = 'coefficient_above_one'
    elif dbin <= tol:
        label = 'numerically_fair'
    else:
        label = 'numerically_nonnegative_nonbinary'
    return {
        'task': task['id'], 'mask': mask, 'bits_by_group': bits,
        'classification': label,
        'min_coefficient': scientific(lo), 'max_coefficient': scientific(hi),
        'max_imaginary_residual': scientific(imag),
        'max_product_target_error': scientific(product_error),
        'max_distance_to_binary': scientific(dbin),
        'rounding_error': scientific(rounding_error),
        'rounded_binary_product_exact': exact_round_product,
        'rounded_A': rounded_a if exact_round_product else None,
        'rounded_B': rounded_b if exact_round_product else None,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--diagnostics', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--tasks', default='18,19,20,21,22,23')
    parser.add_argument('--dps', type=int, default=80)
    parser.add_argument('--tol-exp', type=int, default=55)
    parser.add_argument('--max-steps', type=int, default=80)
    parser.add_argument('--limit-per-task', type=int)
    args = parser.parse_args()
    if args.dps < 60 or not (20 <= args.tol_exp < args.dps-10):
        parser.error('require dps>=60 and 20<=tol-exp<dps-10')
    mp.mp.dps = args.dps
    task_ids = {int(x) for x in args.tasks.split(',')}
    source = json.loads(args.input.read_text())
    prior = json.loads(args.diagnostics.read_text())
    tasks = {x['id']: x for x in source['tasks'] if x['id'] in task_ids}
    chosen = {t: [] for t in task_ids}
    for row in prior:
        t = row.get('task')
        if t in chosen and row.get('classification') == 'nearunknown':
            chosen[t].append(int(row['mask']))
    if args.limit_per_task is not None:
        chosen = {t: v[:args.limit_per_task] for t, v in chosen.items()}
    if set(tasks) != task_ids or any(not chosen[t] for t in task_ids):
        raise ValueError('requested task missing or has no selected nearunknown masks')
    results, root_summary = [], {}
    for t in sorted(task_ids):
        root_data = diagnose_task(tasks[t], chosen[t], args.dps, args.tol_exp, args.max_steps)
        root_summary[str(t)] = {
            k: scientific(v) if isinstance(v, mp.mpf) else v
            for k, v in root_data.items() if k not in ('roots', 'factors', 'root_logs')
        }
        for mask in chosen[t]:
            results.append(diagnose_mask(root_data, tasks[t], mask, root_data.get('tolerance', mp.mpf(0))))
    counts = Counter(x['classification'] for x in results)
    report = {
        'status': 'completed_high_precision_numerical_diagnostic',
        'numeric_only': True,
        'certificate_status': 'not_a_root_allocation_certificate',
        'input': str(args.input), 'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
        'diagnostics': str(args.diagnostics), 'diagnostics_sha256': hashlib.sha256(args.diagnostics.read_bytes()).hexdigest(),
        'tasks': sorted(task_ids), 'dps': args.dps, 'tolerance_exponent': args.tol_exp,
        'max_newton_steps': args.max_steps, 'selected_masks_by_task': {str(k): v for k, v in chosen.items()},
        'root_summary': root_summary, 'classification_counts': dict(counts), 'results': results,
        'mpmath': mp.__version__, 'python': sys.version,
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'masks': len(results), 'counts': dict(counts), 'out': str(args.out)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
