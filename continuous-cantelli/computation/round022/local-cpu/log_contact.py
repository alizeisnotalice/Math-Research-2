"""Adaptive quadrature for the logarithmic approach to asymmetric contact."""
from concurrent.futures import ProcessPoolExecutor
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import dawsn
from scipy.stats import qmc

LS = (8, 12, 18, 28, 48, 80, 128, 224, 384, 640, 1024, 2048)
CROSSCHECK_LS = (18, 128, 1024)
SHARDS = 12
SOBOL_POWER = 10
SEED = 201119
RADIUS = .01
MIN_S = -36
MAX_S = 64
PANEL_NODES = (-18, -10, -5, -2, 0, 2, 4, 7, 11, 18, 28, 42, 52)


def adaptive_side(x, coefficient, power, lam):
    upper = min(MAX_S, lam + math.log(RADIUS))

    def integrand(s):
        t2 = math.exp(2 * s)
        u = math.exp(s - lam)
        var = 1 + coefficient * t2 * (lam - s)**power
        return t2 * math.exp(-t2 / (2 * var) - x * u - u*u/2) / (var * math.sqrt(var))

    return quad(integrand, MIN_S, upper,
                points=[p for p in PANEL_NODES if p < upper],
                epsabs=2e-13, epsrel=2e-13, limit=100)


def mp_side(x, coefficient, power, lam):
    # Independent tanh–sinh adaptive panels, entirely in mp arithmetic.
    x, coefficient, power = (mp.mpf(value) for value in (x, coefficient, power))
    upper = min(MAX_S, lam + mp.log(mp.mpf('0.01')))

    def integrand(s):
        t2 = mp.exp(2 * s)
        u = mp.exp(s - lam)
        var = 1 + coefficient * t2 * (lam - s)**power
        return t2 * mp.exp(-t2 / (2 * var) - x * u - u*u/2) / (var * mp.sqrt(var))

    points = [-36] + [p for p in (-20, -8, -2, 2, 7, 16, 32, 48) if p < upper] + [upper]
    return mp.quad(integrand, points)


def frozen(A):
    z = math.sqrt(0.5 / A)
    return float(2 * z * dawsn(z))


def params():
    raw = qmc.Sobol(d=5, scramble=True, seed=SEED).random_base2(SOBOL_POWER)
    for j, (x_in, c_in, d_in, p_in, q_in) in enumerate(raw):
        x = float(8*x_in - 4)
        c = float(math.exp((2*c_in - 1)*math.log(6)))
        d = float(math.exp((2*d_in - 1)*math.log(6)))
        p = float(0.2 + 3*p_in)
        q = float(0.2 + 3*q_in)
        family = ('tuned', 'equal-power', 'free')[j % 3]
        if family == 'tuned':
            d = 1 / c
        if family != 'free':
            q = p
        yield j, family, x, c, d, p, q


def worker(batch):
    mp.mp.dps = 26
    cpu_start = time.process_time()
    cases = []
    for index, family, x, c, d, p, q in batch:
        rows = []
        for lam in LS:
            plus, plus_err = adaptive_side(x, d, -q, lam)
            minus, minus_err = adaptive_side(-x, c, p, lam)
            scaled = plus - minus
            a_plus, a_minus = d * lam**(-q), c * lam**p
            model = frozen(a_plus) - frozen(a_minus)
            first = a_plus - 1 / a_minus
            cross_error = None
            if lam in CROSSCHECK_LS:
                mp_plus = mp_side(x, d, -q, lam)
                mp_minus = mp_side(-x, c, p, lam)
                cross_error = float(max(abs(mp_plus - plus), abs(mp_minus - minus)))
                if cross_error > 2e-11:
                    raise ArithmeticError((index, lam, cross_error))
            rows.append({
                'log_inverse_h': lam, 'scaled': scaled, 'error_to_limit': scaled - 1,
                'frozen': model, 'frozen_residual': scaled - model,
                'first_prediction': first,
                'quad_error_estimate': plus_err + minus_err, 'cross_error': cross_error,
            })
        cases.append({'index': index, 'family': family, 'x': x, 'c': c, 'd': d,
                      'p': p, 'q': q, 'rows': rows})
    return {'pid': os.getpid(), 'process_cpu_seconds': time.process_time() - cpu_start,
            'cases': cases}


def ps_snapshot(parent_pid):
    process_lines = subprocess.check_output(
        ['/bin/ps', '-axo', 'pid=,ppid=,%cpu=,time=,state='], text=True
    ).splitlines()
    return [row for line in process_lines if (row := line.split())[1] == str(parent_pid)]


def main():
    configurations = list(params())
    samplings = []
    start = time.monotonic()
    pid = os.getpid()
    snapshot_re = re.compile(r'CPU usage: (.*?)% user, (.*?)% sys, (.*?)% idle')
    with ProcessPoolExecutor(max_workers=SHARDS) as pool:
        pending = [pool.submit(worker, configurations[i::SHARDS]) for i in range(SHARDS)]
        print(f'parent={pid}, workers={SHARDS}, cases={len(configurations)}', flush=True)
        while not all(future.done() for future in pending):
            report = subprocess.check_output(
                ['/usr/bin/top', '-l', '2', '-s', '1', '-n', '0'], text=True
            )
            whole = snapshot_re.findall(report)[-1]
            sample = {
                'elapsed_s': time.monotonic() - start,
                'whole_cpu_percent_user_sys_idle': list(map(float, whole)),
                'ps_pid_ppid_percent_cputime_state': ps_snapshot(pid),
            }
            samplings.append(sample)
            print(json.dumps(sample), flush=True)
            time.sleep(.2)
        shards = [future.result() for future in pending]
    worker_pids = [shard['pid'] for shard in shards]
    assert len(set(worker_pids)) == SHARDS
    active = set(map(int, subprocess.check_output(['/bin/ps', '-axo', 'pid='], text=True).split()))
    assert not active.intersection(worker_pids)
    result = {
        'model': 'v(x-u)/u^2=c log(1/u)^p, v(x+u)/u^2=d log(1/u)^(-q)',
        'normalization': 'h D1 / (p1(x)/sqrt(2*pi)), limit 1 for all tested cases',
        'radius': RADIUS, 'log_t_clip': [MIN_S, MAX_S], 'quad_epsabs_epsrel': 2e-13,
        'mpmath_dps': 26, 'qmc_seed': SEED, 'python': sys.version,
        'compiler': platform.python_compiler(), 'platform': platform.platform(),
        'numpy': np.__version__, 'scipy': scipy.__version__, 'mpmath': mp.__version__,
        'elapsed_seconds': time.monotonic() - start, 'workers_verified_exited': True,
        'samplings': samplings, 'shards': shards,
    }
    Path('log-contact-results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(f"done. {len(configurations)} cases; workers verified exited, {result['elapsed_seconds']:.1f}s.",
          flush=True)


if __name__ == '__main__':
    main()
