"""Higher signed Peano perturbations about asymmetric quartic contact."""
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
SEED = 430118
H_EXPONENTS = (7, 10, 13, 17, 20)
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
    samples = qmc.Sobol(d=6, scramble=True, seed=SEED).random_base2((n-1).bit_length())[:n]
    for j, (x, c_val, plus_val, minus_val, h_in, b_in) in enumerate(samples):
        beta = (8.5, 9.0, 10.0, 11.5, 12.5)[j % 5]
        c = .35 * math.exp(float(c_val) * math.log(3/.35))
        plus = float(1.05*plus_val-.525)
        minus = float(1.05*minus_val-.525)
        yield j, beta, float(6*x-3), c, plus, minus, float(1.8*h_in-.9), float(b_in-.5)


def cuts(top):
    return [p for p in LOG_PANELS if p < top] + [top]


def high_precision(case, k, exponent):
    _, beta, x, c, plus, minus, H1, H2 = (mp.mpf(entry) for entry in case)
    h = mp.mpf(2)**(-exponent)
    end = mp.log(mp.mpf(RADIUS)/h)
    scale = c/mp.mpf(RADIUS)**(beta-4)

    def integrand(log_t):
        t = mp.exp(log_t)
        u = h*t
        common = t * mp.exp(-u*u/2)
        change = scale*u**beta/h**2
        pieces = []
        for sigma, d in ((1,plus),(-1,minus)):
            ref = 1+c*u**4*(1+sigma*H1*u+H2*u*u)**2/h**2
            shifted = kernel(k,t,ref+d*change,mp.exp,mp.sqrt)
            base = kernel(k,t,ref,mp.exp,mp.sqrt)
            pieces.append(sigma**k*mp.exp(-sigma*x*u)*(shifted-base))
        return common*mp.fsum(pieces)

    return h**(2-beta)*mp.quad(integrand,cuts(end))/scale


def scipy_perturbation(case, k, exponent):
    _, beta, x, c, plus, minus, H1, H2 = case
    h = 2.0**(-exponent)
    end = math.log(RADIUS/h)
    scale = c/RADIUS**(beta-4)

    def side(t, u, sigma, d):
        v0 = 1 + c*u**4*(1+sigma*H1*u+H2*u*u)**2/h**2
        dv = scale*d*u**beta/h**2
        variances = (v0 + (q+1)*dv/2 for q in G_NODES)
        j = math.fsum(w*kernel(k+2,t,var,math.exp,math.sqrt)
                      for w,var in zip(G_WEIGHTS,variances))
        return d*j/4

    def integrand(log_t):
        t = math.exp(log_t)
        u = h*t
        return t**(beta+1)*math.exp(-u*u/2)*(
            math.exp(-x*u)*side(t,u,1,plus)
            +(-1)**k*math.exp(x*u)*side(t,u,-1,minus)
        )

    panels = cuts(end)
    return quad(integrand,panels[0],panels[-1],points=panels[1:-1],
                epsabs=1e-8,epsrel=1e-12,limit=120)


def coefficient(case, k):
    _, beta, _, _, plus, minus, _, _ = case
    parity = k % 2
    moment = 2**((beta+parity-1)/2)*math.gamma((beta+parity+1)/2)
    factors = math.prod(beta-j for j in range(parity, k+1, 2))
    return (plus + (-1)**k*minus)*moment*factors/2


def worker(batch, smoke):
    mp.mp.dps = 92
    started = utc()
    cpu_start = time.process_time()
    records = []
    for case in batch:
        entry = dict(zip(("index", "beta", "x", "c", "d_plus", "d_minus", "H1", "H2"), case))
        series = []
        for k in (int(case[1]), int(case[1])+1):
            expected = coefficient(case, k)
            for exponent in ((7, 20) if smoke else H_EXPONENTS):
                mp_value = high_precision(case, k, exponent)
                fast, estimate = scipy_perturbation(case, k, exponent)
                error = abs(float(mp_value)-fast)
                if error > 5e-8*max(1, abs(fast)):
                    raise ArithmeticError((case, k, exponent, float(mp_value), fast, error))
                row = {
                    "k": k, "h_two_exponent": -exponent, "coefficient": expected,
                    "scaled_difference": mp.nstr(mp_value, 36),
                    "coefficient_residual": float(mp_value)-expected,
                    "scipy_cross_error": error, "scipy_error_estimate": estimate,
                }
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
    parser.add_argument("--cases", type=int, default=48)
    parser.add_argument("--out", default="peano-results.json")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--start-utc")
    args = parser.parse_args()
    rows = list(parameters(args.cases))
    if args.smoke:
        rows = rows[:2]
    workers = 2 if args.smoke else WORKERS
    if args.start_utc:
        target = datetime.fromisoformat(args.start_utc)
        while (remain := (target - datetime.now(timezone.utc)).total_seconds()) > 0:
            time.sleep(min(remain, .2))
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
        "radius": RADIUS, "seed": SEED, "model": "v0=c*u^4*(1+H1*u+H2*u^2)^2; v=v0+scale*d_side*abs(u)^beta",
        "normalization": "difference times h^(k+2-beta)/(scale*p1(x)/sqrt(2*pi)), scale=c/r^(beta-4)",
        "versions": {"python": sys.version, "platform": platform.platform(),
                     "numpy": np.__version__, "mpmath": mp.__version__, "scipy": scipy.__version__},
        "main_dps": 92, "worker_ps_exit_verified": True,
        "sampling": samples, "shards": results,
    }
    Path(args.out).write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"finished_utc": finished_utc, "elapsed_s": result["elapsed_s"],
                      "workers_verified_exited": True}), flush=True)


if __name__ == "__main__":
    main()
