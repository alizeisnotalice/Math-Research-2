#!/usr/bin/env python3
"""Exact scalar checks for the analytic finite all-prefix area certificate."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.perf_counter()
    spec = json.loads(args.input.read_text())
    n, m = spec['n'], F(spec['total_mass'])
    intervals = [(F(a), F(b)) for a, b in spec['axis_intervals']]
    assert n == 400001 and n % 2 == 1 and m == F(99, 100)
    assert intervals == [(F(-49, 100), F(49, 100)), (F(51, 100), F(149, 100))]
    assert spec['rectangle_error_constant'] == 8 and spec['prefix_error_constant'] == 16
    # All bounds on delta below are strict on the open intervals.
    delta_lo, delta_hi = F(1, 100), F(99, 100)
    sqrt_m_lower = F(99, 100)
    r_lower = sqrt_m_lower / 2
    assert 0 < sqrt_m_lower and sqrt_m_lower**2 < m < 1
    assert r_lower > F(49, 100)
    factor_lower = r_lower + F(1, 2) - delta_hi
    factor_upper = F(1, 2) + F(1, 2) - delta_lo
    a_lower = factor_lower**2 - F(8, n)
    a_upper = factor_upper**2 + F(8, n)
    gamma_lower = sqrt_m_lower - 1 + 2 * delta_lo
    active_distance_lower = ((1 - m) - F(16, n)) / 2
    cap_slack_lower = active_distance_lower * gamma_lower - F(16, n)
    assert factor_lower == F(1, 200) and factor_upper == F(99, 100)
    assert a_lower > 0 and a_upper < m
    assert gamma_lower == F(1, 100) and active_distance_lower > 0
    assert cap_slack_lower == F(1, 20000) - F(402, 25 * n) > 0
    width = 1 - F(1, n*n)
    left_margin = min(a for a, _ in intervals) + F(1, 2)
    right_margin = width + F(1, 2) - max(b for _, b in intervals)
    axis_length = sum(b-a for a, b in intervals)
    area = axis_length**2
    assert left_margin > 0 and right_margin > 0
    assert area == F(2401, 625) > F(7, 2)
    result = dict(
        status='PASS_EXACT_SCALAR_FINITE_ALL_PREFIX_AREA_CERTIFICATE',
        source_sha256=sha(Path(__file__)), input_sha256=sha(args.input),
        python=sys.version, platform=platform.platform(), n=n,
        formula_specified_old_sources=n*n, total_mass=str(m),
        radius='sqrt(99/100)/2',
        strict_bounds={key: str(value) for key, value in dict(
            sqrt_mass_lower=sqrt_m_lower, radius_lower=r_lower,
            factor_lower=factor_lower, factor_upper=factor_upper,
            capture_mass_lower=a_lower, capture_mass_upper=a_upper,
            gamma_lower=gamma_lower, active_distance_lower=active_distance_lower,
            all_active_prefix_cap_slack_lower=cap_slack_lower,
            Cn_left_margin=left_margin, Cn_right_margin=right_margin).items()},
        area_lower_certificate=str(area), area_above_7_over_2=str(area-F(7, 2)),
        elapsed_seconds=time.perf_counter()-start,
        scope='Fraction scalar verification plus section-5 analytic uniform rectangle error; no source/subset enumeration, no actual halo/W/Wprime/deletion calculation, no new C2 bound.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
