"""Cross-check the fourth derivative for the exact local linear–cubic model."""
from concurrent.futures import ProcessPoolExecutor
from itertools import pairwise
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import mpmath as mp

# Rows give (a, d), where a = 1/c and x = 3d - 2a.
CASES = (
    ("unit", "1", "1"),
    ("left-heavy", "4", "0.2"),
    ("near-negative", "1", "0.125"),
    ("near-positive", "1", "0.25"),
    ("origin", "0.75", "0.5"),
    ("sharper-cubic", "0.125", "4"),
    ("narrow-linear", "6", "0.8"),
    ("wide-linear", "0.125", "0.125"),
    ("small-positive", "0.05", "0.075"),
    ("root-a=d", "root", "root"),
    ("negative-tail", "2", "0.04"),
    ("positive-x", "0.3", "0.75"),
)
EXPONENTS = (12, 15, 18, 22, 26)


def case_parameters(case):
    label, a_string, d_string = case
    if a_string == "root":
        a = d = mp.sqrt(mp.mpf(27) / 6955)
    else:
        a, d = mp.mpf(a_string), mp.mpf(d_string)
    return label, a, d, 3 * d - 2 * a


def geometric_cuts(last, base, start=0):
    cuts = [mp.mpf(start)]
    current = start or 1
    while current < last:
        if current > cuts[-1]:
            cuts.append(mp.mpf(current))
        current *= base
    return cuts + [last]


def k4(z, variance):
    scaled_sq = z * z / variance
    return (3 - 6 * scaled_sq + scaled_sq**2) * mp.exp(-scaled_sq / 2) / (
        variance**2 * mp.sqrt(variance)
    )


def both_sides(a, d, x, h, method, dps):
    with mp.workdps(dps):
        r = mp.mpf("0.005") / max(d, 1)
        left_cut = r / (a * h**2)
        right_cut = r / h
        left_points = geometric_cuts(left_cut, 4)

        core_points = [mp.mpf(j) for j in range(0, 41, 4) if j < right_cut]
        if right_cut > 40:
            right_points = core_points[:-1] + geometric_cuts(right_cut, 4, 40)
        else:
            right_points = core_points + [right_cut]

        def left(t):
            u = a * h**2 * t
            return mp.exp(x * u - u**2 / 2) * k4(a * h * t, 1 + t)

        def right(t):
            u = h * t
            return mp.exp(-x * u - u**2 / 2) * k4(t, 1 + d * h * t**3)

        def integrate(function, points):
            return mp.fsum(mp.quad(function, [lo, hi], method=method)
                           for lo, hi in pairwise(points))

        return a * integrate(left, left_points) / h**2, integrate(right, right_points) / h**3


def worker(case):
    mp.mp.dps = 75
    label, a, d, x = case_parameters(case)
    s_left = 4 * x * a**2 - 8 * a**3
    s_right = 3 * x - x**3 + 30 * d * (x**2 - 1) - 630 * d**2 * x + 7560 * d**3
    expected = s_left + s_right
    rows = []
    cpu_start = time.process_time()
    for k in EXPONENTS:
        h = mp.mpf(2) ** -k
        left, right = both_sides(a, d, x, h, "gauss-legendre", 60)
        reference_left, reference_right = both_sides(a, d, x, h, "tanh-sinh", 75)
        total = left + right
        cross_error = abs(total - reference_left - reference_right)
        if cross_error > mp.mpf("1e-24") * max(1, abs(total)):
            raise ArithmeticError((label, k, mp.nstr(cross_error, 24)))
        rows.append({
            "h_exponent": k,
            "epsilon_exponent": 2 * k,
            "total": mp.nstr(total, 36),
            "error_to_coefficient": mp.nstr(total - expected, 30),
            "left_after_leading_removed": mp.nstr(left - 2 * a / h**2, 30),
            "right_after_leading_removed": mp.nstr(right - (x - 3 * d) / h**2, 30),
            "cross_quadrature_absolute_difference": mp.nstr(cross_error, 20),
        })
    return {
        "label": label,
        "pid": os.getpid(),
        "cpu_seconds": time.process_time() - cpu_start,
        "a": mp.nstr(a, 36), "c": mp.nstr(1 / a, 36),
        "d": mp.nstr(d, 36), "x": mp.nstr(x, 36),
        "left_expected": mp.nstr(s_left, 36),
        "right_expected": mp.nstr(s_right, 36),
        "predicted_S": mp.nstr(expected, 36),
        "rows": rows,
    }


def main():
    samples = []
    start = time.monotonic()
    with ProcessPoolExecutor(max_workers=12) as pool:
        pending = [pool.submit(worker, case) for case in CASES]
        print(f"parent={os.getpid()} workers=12", flush=True)
        while not all(future.done() for future in pending):
            proc_table = subprocess.check_output(
                ["/bin/ps", "-axo", "pid=,ppid=,%cpu=,time=,state="], text=True
            )
            rows = [line.split() for line in proc_table.splitlines()]
            children = [row for row in rows if int(row[1]) == os.getpid()]
            samples.append({
                "seconds": round(time.monotonic() - start, 3),
                "columns": ["pid", "ppid", "cpu_percent", "cpu_time", "state"],
                "children": children,
            })
            time.sleep(1)
        cases = [future.result() for future in pending]

    finished = set(str(case["pid"]) for case in cases)
    active = set(subprocess.check_output(["/bin/ps", "-axo", "pid="], text=True).split())
    assert len(finished) == 12 and not finished.intersection(active)
    output = {
        "model": "v(x-u)=c*u, v(x+u)=d*u^3, 0<u<r, x=3d-2/c",
        "normalization": "h*D4/(p1(x)/sqrt(2*pi)), h=sqrt(epsilon)",
        "cutoff": "r=0.005/max(1,d)",
        "methods": {"gauss-legendre_dps": 60, "tanh-sinh_dps": 75},
        "relative_or_absolute_agreement_threshold": "1e-24",
        "python": sys.version,
        "python_compiler": platform.python_compiler(),
        "mpmath": mp.__version__,
        "platform": platform.platform(),
        "logical_cores": os.cpu_count(),
        "worker_count": 12,
        "elapsed_seconds": time.monotonic() - start,
        "worker_exit_ps_verified": True,
        "ps_samples": samples,
        "cases": cases,
    }
    Path("quadrature-results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(f"done. {len(cases)} cases, {len(samples)} ps samples, workers verified exited",
          flush=True)


if __name__ == "__main__":
    main()
