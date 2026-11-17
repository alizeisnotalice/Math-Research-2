#!/usr/bin/env python3
"""Independent CPU reconstruction verifier for native mask-search candidates.

It reads the original root-group task, never LR/LI, DC/DS, GPU scores, Swift
code, or a production success flag.  Root data are floating point: this is a
numerical candidate diagnostic, never an algebraic certificate.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any


DEFAULT_SOURCE = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-macbook-gpu52-direct/inputs.json')
DEFAULT_MASKS = (0, 1, 17, 1234, 23166975, 23166976)


def finite_complex(z: complex) -> bool:
    return math.isfinite(z.real) and math.isfinite(z.imag)


def convolve(left: list[complex], right: list[complex]) -> list[complex]:
    out = [0j] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            out[i + j] += x * y
    return out


def group_factor(group: Any) -> list[complex]:
    if not isinstance(group, list) or not group:
        raise ValueError('each root group must be a nonempty list')
    coeff = [1 + 0j]
    for pair in group:
        if not (isinstance(pair, list) and len(pair) == 2):
            raise ValueError('a root must be [real,imag]')
        z = complex(float(pair[0]), float(pair[1]))
        if not finite_complex(z) or z == 0:
            raise FloatingPointError('root is nonfinite or zero')
        # prod (x-z)/(-z) = prod (1-x/z), so the constant coefficient is one.
        coeff = convolve(coeff, [1 + 0j, -1 / z])
    return coeff


def task_from_source(doc: dict[str, Any], task_number: int) -> dict[str, Any]:
    tasks = doc.get('tasks')
    if not isinstance(tasks, list) or not (0 <= task_number < len(tasks)):
        raise ValueError('task index outside source task list')
    task = tasks[task_number]
    if not isinstance(task, dict) or task.get('id') != task_number:
        raise ValueError('task id/index mismatch')
    groups = task.get('groups')
    roots, target = task.get('roots'), task.get('C')
    if not isinstance(groups, int) or groups < 2 or groups >= 32:
        raise ValueError('groups must be an integer in [2,31]')
    if not isinstance(roots, list) or len(roots) != groups:
        raise ValueError('roots group dimension mismatch')
    if not isinstance(target, list) or not target:
        raise ValueError('target C is absent')
    if any(not isinstance(x, (int, float)) or not math.isfinite(float(x)) for x in target):
        raise FloatingPointError('target has nonfinite value')
    if any(float(x) not in (0.0, 1.0) for x in target):
        raise ValueError('target C must be binary')
    return task


def bits_for_mask(groups: int, mask: int) -> list[int]:
    if not isinstance(mask, int) or mask < 0 or mask >= (1 << (groups - 1)):
        raise ValueError(f'mask {mask!r} outside [0,2^(groups-1))')
    # Group zero is the fixed control group in B.  This order matches the
    # documented fixed-width base-two convention but is implemented here anew.
    return [0] + [int((mask >> (groups - 1 - g)) & 1) for g in range(1, groups)]


def serialise(values: list[complex]) -> dict[str, list[float]]:
    return {
        'real': [float(z.real) for z in values],
        'imag': [float(z.imag) for z in values],
    }


def reconstruct(task: dict[str, Any], mask: int, real_tol: float, boundary_tol: float, product_tol: float) -> dict[str, Any]:
    groups = task['groups']
    try:
        bits = bits_for_mask(groups, mask)
        factors = [group_factor(group) for group in task['roots']]
        A, B = [1 + 0j], [1 + 0j]
        for bit, factor in zip(bits, factors):
            if bit:
                A = convolve(A, factor)
            else:
                B = convolve(B, factor)
        AB = convolve(A, B)
        target = [complex(float(x), 0.0) for x in task['C']]
        if len(AB) != len(target):
            raise ValueError(f'product degree dimension {len(AB)} != target {len(target)}')
        all_values = A + B + AB
        if not all(finite_complex(z) for z in all_values):
            return {'mask': mask, 'classification': 'nonfinite', 'reason': 'nonfinite reconstructed coefficient'}
        imag_residual = max(abs(z.imag) for z in all_values)
        product_error = max(abs(z - t) for z, t in zip(AB, target))
        ar, br = [z.real for z in A], [z.real for z in B]
        min_coeff, max_coeff = min(ar + br), max(ar + br)
        distance_to_binary = max(min(abs(v), abs(v - 1.0)) for v in ar + br)
        # The source roots are decimal approximations.  The all-groups product
        # residual against C calibrates, but does not certify, an error scale.
        effective_boundary_tol = max(boundary_tol, 2.0 * product_error)
        numeric_score = max(0.0, -min_coeff, max_coeff - 1.0, imag_residual)
        if imag_residual > real_tol:
            classification, reason = 'nearunknown', 'root reconstruction is not sufficiently real'
        elif product_error > product_tol:
            classification, reason = 'nearunknown', 'root reconstruction is not sufficiently target-consistent'
        elif min_coeff < -effective_boundary_tol:
            classification, reason = 'negative', 'a reconstructed factor has a strictly negative coefficient'
        elif max_coeff > 1.0 + effective_boundary_tol:
            classification, reason = 'nonuniform', 'a reconstructed nonnegative candidate exceeds coefficient one'
        elif distance_to_binary <= effective_boundary_tol:
            classification, reason = 'fair', 'both reconstructed factors are numerically binary at the calibrated float scale'
        elif min_coeff >= -effective_boundary_tol and max_coeff <= 1.0 + effective_boundary_tol:
            classification, reason = 'nonuniform', 'coefficients are numerically within [0,1] but not binary'
        else:
            classification, reason = 'nearunknown', 'boundary comparison is numerically indeterminate'
        return {
            'mask': mask,
            'bits_by_group': bits,
            'classification': classification,
            'reason': reason,
            'finite': True,
            'degree_A': len(A) - 1,
            'degree_B': len(B) - 1,
            'min_coefficient': min_coeff,
            'max_coefficient': max_coeff,
            'max_imaginary_residual': imag_residual,
            'max_product_target_error': product_error,
            'max_distance_to_binary': distance_to_binary,
            'effective_boundary_tolerance': effective_boundary_tol,
            'numeric_boundary_score': numeric_score,
            'A': serialise(A),
            'B': serialise(B),
            'AB': serialise(AB),
        }
    except FloatingPointError as exc:
        return {'mask': mask, 'classification': 'nonfinite', 'reason': str(exc), 'finite': False}
    except (ValueError, TypeError, OverflowError) as exc:
        return {'mask': mask, 'classification': 'invalid', 'reason': str(exc), 'finite': False}


def run_self_tests(source: Path) -> dict[str, Any]:
    doc = json.loads(source.read_text())
    task = task_from_source(doc, 2)
    # The saved roots are decimal approximations; the product residual sets
    # the adaptive float scale for the all-B control.
    strict_control = reconstruct(task, 0, 1e-8, 1e-8, 1e-3)
    assert strict_control['classification'] == 'fair', strict_control
    synthetic_fair = {
        'groups': 2,
        'roots': [[[-1.0, 0.0]], [[0.0, 1.0], [0.0, -1.0]]],
        'C': [1, 1, 1, 1],
    }
    fair = reconstruct(synthetic_fair, 0, 1e-12, 1e-12, 1e-12)
    assert fair['classification'] == 'fair', fair
    bad_mask = reconstruct(task, 1 << (task['groups'] - 1), 1e-8, 1e-8, 1e-7)
    assert bad_mask['classification'] == 'invalid', bad_mask
    malformed = copy.deepcopy(task)
    malformed['roots'] = malformed['roots'][:-1]
    try:
        # Revalidate through a temporary source-like document to exercise the dimension guard.
        task_from_source({'tasks': [malformed, malformed, malformed]}, 2)
        raise AssertionError('dimension mismatch unexpectedly accepted')
    except ValueError:
        pass
    nan_task = copy.deepcopy(task)
    nan_task['roots'][0][0][0] = float('nan')
    nan_case = reconstruct(nan_task, 1, 1e-8, 1e-8, 1e-7)
    assert nan_case['classification'] == 'nonfinite', nan_case
    return {
        'status': 'passed_independent_negative_tests',
        'decimal_control_mask_at_calibrated_scale': strict_control['classification'],
        'synthetic_exact_fair_mask': fair['classification'],
        'out_of_range_mask': bad_mask['classification'],
        'dimension_mismatch': 'rejected',
        'nan_root': nan_case['classification'],
    }


def parse_masks(text: str) -> list[int]:
    return [int(part) for part in text.split(',') if part.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--task', type=int, default=2)
    parser.add_argument('--masks', default=','.join(map(str, DEFAULT_MASKS)))
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--real-tol', type=float, default=1e-8)
    parser.add_argument('--boundary-tol', type=float, default=1e-8)
    parser.add_argument('--product-tol', type=float, default=1e-3)
    args = parser.parse_args()
    if min(args.real_tol, args.boundary_tol, args.product_tol) <= 0:
        parser.error('tolerances must be positive')
    doc = json.loads(args.source.read_text())
    task = task_from_source(doc, args.task)
    masks = parse_masks(args.masks)
    results = [reconstruct(task, mask, args.real_tol, args.boundary_tol, args.product_tol) for mask in masks]
    summary: dict[str, int] = {}
    for row in results:
        summary[row['classification']] = summary.get(row['classification'], 0) + 1
    report = {
        'status': 'completed_independent_cpu_reconstruction',
        'numeric_only': True,
        'certificate_status': 'not_a_certificate_float_root_input',
        'source': str(args.source),
        'source_sha256': hashlib.sha256(args.source.read_bytes()).hexdigest(),
        'task': args.task,
        'groups': task['groups'],
        'masks': masks,
        'tolerances': {'real': args.real_tol, 'boundary': args.boundary_tol, 'product': args.product_tol},
        'classification_counts': summary,
        'results': results,
        'negative_tests': run_self_tests(args.source) if args.self_test else None,
        'python': sys.version,
        'source_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': report['status'], 'classification_counts': summary, 'out': str(args.out)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
