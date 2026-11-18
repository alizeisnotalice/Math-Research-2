"""Adaptive core integrals for smooth quartic baselines and signed perturbations."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
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
from scipy.stats import qmc

RADIUS = 0.125
SEED = 230118
H_EXPONENTS = (8, 12, 17, 23, 29, 35)
BASELINE_EXPONENTS = (8, 17, 29, 35)
K_ORDERS = (4, 5)
LOG_PANELS = (-54, -32, -19, -11, -6, -3, -1, 0, 1, 2, 3, 4, 6, 9, 14, 20, 27)
WORKERS = 12
CPU_PATTERN = re.compile(r"CPU usage: (.*?)% user, (.*?)% sys, (.*?)% idle")
G_NODES, G_WEIGHTS = (tuple(map(float, arr)) for arr in np.polynomial.legendre.leggauss(12))


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def hermite(order, z):
    previous, current = 1, z
    for j in range(1, order):
        previous, current = current, z * current - j * previous
    return previous if order == 0 else current


def kernel(order, t, var, exp, sqrt):
    z = t / sqrt(var)
    return hermite(order, z) * exp(-z*z/2) / (var**((order+1)//2) * (sqrt(var) if order % 2 == 0 else 1))


def parameters(n):
    samples = qmc.Sobol(d=4, scramble=True, seed=SEED).random_base2(int(math.log2(n)))
    for j, (x, c_val, plus_val, minus_val) in enumerate(samples):
        beta = 4.5 if j % 2 == 0 else 5.0
        c = .35 * math.exp(float(c_val) * math.log(3/.35))
        plus = float(1.05*plus_val - .525) * c / RADIUS**(beta-4)
        minus = float(1.05*minus_val - .525) * c / RADIUS**(beta-4)
        yield j, beta, float(6*x-3), c, plus, minus


def cuts(top):
    return [p for p in LOG_PANELS if p < top] + [top]


def high_precision(case, k, exponent, baseline=False):
    _, beta, x, c, plus, minus = (mp.mpf(entry) for entry in case)
    h = mp.mpf(2)**(-exponent)
    end = mp.log(mp.mpf(RADIUS)/h)
    zero_adjust = h**(-k) if baseline else h**(2-beta)

    def integrand(log_t):
        t = mp.exp(log_t)
        u = h*t
        v0 = 1 + c*u**4/h**2
        base_f = kernel(k, t, v0, mp.exp, mp.sqrt)
        common = t * mp.exp(-u*u/2)
        if baseline:
            # The signed source term does not lose precision at small u.
            both = 2*(mp.cosh(x*u) if k % 2 == 0 else -mp.sinh(x*u))
            return common * both * (base_f - kernel(k, t, 1, mp.exp, mp.sqrt))
        du = u**beta/h**2
        fp = kernel(k, t, v0+plus*du, mp.exp, mp.sqrt)
        fm = kernel(k, t, v0+minus*du, mp.exp, mp.sqrt)
        return common * (mp.exp(-x*u)*(fp-base_f) + (-1)**k*mp.exp(x*u)*(fm-base_f))

    return zero_adjust*mp.quad(integrand, cuts(end))


def scipy_perturbation(case, k, exponent):
    _, beta, x, c, plus, minus = case
    h = 2.0**(-exponent)
    end = math.log(RADIUS/h)

    def side(t, u, v0, d):
        dv = d*u**beta/h**2
        # Integrate the heat-equation difference quotient over its variance path.
        variances = (v0 + (q+1)*dv/2 for q in G_NODES)
        j = math.fsum(w*kernel(k+2, t, var, math.exp, math.sqrt)
                      for w, var in zip(G_WEIGHTS, variances))
        return d*j/4

    def integrand(log_t):
        t = math.exp(log_t)
        u = h*t
        v0 = 1 + c*u**4/h**2
        return t**(beta+1)*math.exp(-u*u/2) * (
            math.exp(-x*u)*side(t, u, v0, plus) +
            (-1)**k*math.exp(x*u)*side(t, u, v0, minus)
        )

    panels = cuts(end)
    return quad(integrand, panels[0], panels[-1], points=panels[1:-1],
                epsabs=1e-10, epsrel=1e-12, limit=120)


def coefficient(case, k):
    _, beta, _, _, plus, minus = case
    parity = k % 2
    moment = 2**((beta+parity-1)/2)*math.gamma((beta+parity+1)/2)
    factors = math.prod(beta-j for j in range(parity, k+1, 2))
    return (plus + (-1)**k*minus)*moment*factors/2


def worker(batch, smoke):
    mp.mp.dps = 77
    started = utc()
    cpu_start = time.process_time()
    records = []
    for case in batch:
        entry = dict(zip(("index", "beta", "x", "c", "d_plus", "d_minus"), case))
        series = []
        for k in K_ORDERS:
            expected = coefficient(case, k)
            for exponent in ((8, 35) if smoke else H_EXPONENTS):
                mp_value = high_precision(case, k, exponent)
                fast, estimate = scipy_perturbation(case, k, exponent)
                error = abs(float(mp_value)-fast)
                if error > 1e-8*max(1, abs(fast)):
                    raise ArithmeticError((case, k, exponent, float(mp_value), fast, error))
                row = {
                    "k": k, "h_two_exponent": -exponent, "coefficient": expected,
                    "scaled_difference": mp.nstr(mp_value, 36),
                    "coefficient_residual": float(mp_value)-expected,
                    "scipy_cross_error": error, "scipy_error_estimate": estimate,
                }
                if exponent in BASELINE_EXPONENTS:
                    m = high_precision(case, k, exponent, baseline=True)
                    row["unnormalized_scale_baseline"] = mp.nstr(m, 34)
                    if entry["index"] % 32 == 0 and exponent == 35:
                        with mp.workdps(108):
                            more = high_precision(case, k, exponent, baseline=True)
                        baseline_error = abs(m-more)
                        if baseline_error > mp.mpf("1e-14") * max(1, abs(more)):
                            raise ArithmeticError((case, k, exponent, str(baseline_error)))
                        row["baseline_high_dps_error"] = mp.nstr(baseline_error, 22)
                series.append(row)
        entry["series"] = series
        records.append(entry)
    finished = utc()
    return {"pid": os.getpid(), "started_utc": started, "finished_utc": finished,
            "cpu_seconds": time.process_time()-cpu_start, "records": records}


def process_table(pid):
    lines = subprocess.check_output(
        ["/bin/ps", "-axo", "pid=,ppid=,%cpu=,time=,state="], text=True
    ).splitlines()
    return [parts for line in lines if (parts := line.split())[1] == str(pid)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=256)
    parser.add_argument("--out", default="quartic-results.json")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    rows = list(parameters(args.cases))
    if args.smoke:
        rows = rows[:2]
    workers = 2 if args.smoke else WORKERS
    pid = os.getpid()
    start = time.monotonic()
    start_utc = utc()
    print(json.dumps({"parent": pid, "started_utc": start_utc, "cases": len(rows),
                      "workers": workers}), flush=True)
    samples = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(worker, rows[i::workers], args.smoke) for i in range(workers)]
        while not all(f.done() for f in futures):
            top = subprocess.check_output(
                ["/usr/bin/top", "-l", "2", "-s", "1", "-n", "0"], text=True)
            totals = CPU_PATTERN.findall(top)
            sample = {"utc": utc(), "elapsed_s": time.monotonic()-start,
                      "cpu_pct_user_sys_idle": list(map(float, totals[-1])),
                      "ps_pid_ppid_cpu_pct_time_state": process_table(pid)}
            samples.append(sample)
            print(json.dumps(sample), flush=True)
            time.sleep(.15)
        results = [f.result() for f in futures]
    finished_utc = utc()
    worker_ids = {result["pid"] for result in results}
    active = set(map(int, subprocess.check_output(["/bin/ps", "-axo", "pid="], text=True).split()))
    assert len(worker_ids) == workers and not active.intersection(worker_ids)
    result = {
        "start_utc": start_utc, "finish_utc": finished_utc, "elapsed_s": time.monotonic()-start,
        "radius": RADIUS, "seed": SEED, "model": "v0=c*u^4; v=v0+d_side*u^beta",
        "normalization": "difference times h^(k+2-beta)/(p1(x)/sqrt(2*pi)) for perturbation; no h multiplier for baseline",
        "versions": {"python": sys.version, "platform": platform.platform(),
                     "numpy": np.__version__, "mpmath": mp.__version__, "scipy": scipy.__version__},
        "main_dps": 77, "baseline_crosscheck_dps": 108, "worker_ps_exit_verified": True,
        "sampling": samples, "shards": results,
    }
    Path(args.out).write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"finished_utc": finished_utc, "elapsed_s": result["elapsed_s"],
                      "workers_verified_exited": True}), flush=True)


if __name__ == "__main__":
    main()
