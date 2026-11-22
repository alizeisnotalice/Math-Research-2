"""Exact checks of the round-4 smooth cutoff. Uses all available logical CPUs."""
import argparse
from fractions import Fraction
import os
from pathlib import Path
import platform
import subprocess
import sys
from time import monotonic


def scaled_cutoff(x, limit):
    """Return 4*limit times the smooth cutoff and its first derivative."""
    if x <= limit:
        return x**4, 4*x**3
    excess = x - limit
    return (
        limit**4 + 4*limit**3*excess + 6*limit**2*excess**2 + 4*limit*excess**3,
        4*limit**3 + 12*limit**2*excess + 12*limit*excess**2,
    )


def grid(worker, workers):
    checked = 0
    for limit in range(worker + 1, 192 + 1, workers):
        values = [scaled_cutoff(x, limit) for x in range(769)]
        for x in range(385):
            fx, dx = values[x]
            for y in range(x, 385):
                fy, dy = values[y]
                s = x + y
                fs, ds = values[s]
                shifted, _ = values[max(s-1, 0)]
                assert 0 <= fs - fx - fy - x*dy - y*dx <= 6*x*x*y*y
                assert ds >= dx + dy
                assert 0 <= shifted - fs + ds <= 6*s*s
                checked += 1
    return checked


def orbit_checks():
    count = 0
    for p in (
        [Fraction(4, 5), 0, Fraction(1, 5)],
        [Fraction(8, 17), Fraction(1, 2), 0, Fraction(1, 34)],
    ):
        for _ in range(6):
            nxt = [0] * (2*len(p)-2)
            for x, px in enumerate(p):
                for y, py in enumerate(p):
                    nxt[max(x+y-1, 0)] += px*py
            weighted = [2**x*px for x, px in enumerate(p)]
            g = sum(weighted)
            assert sum(x*ax for x, ax in enumerate(weighted)) == g
            t2 = sum(x*x*ax for x, ax in enumerate(weighted))
            for limit in (1, 3, 7, 16, 64, 128, 256):
                before = sum(scaled_cutoff(x, limit)[0]*ax for x, ax in enumerate(weighted))
                after = sum(scaled_cutoff(x, limit)[0]*2**x*px for x, px in enumerate(nxt))
                assert after <= g*before + 3*t2**2 + 12*g*t2
                count += 1
            p = nxt
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=int)
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()
    if args.worker is not None:
        print(grid(args.worker, args.workers))
        return
    start = monotonic()
    script = str(Path(__file__).resolve())
    workers = [
        subprocess.Popen(
            [sys.executable, "-B", script, "--worker", str(i), "--workers", str(args.workers)],
            stdout=subprocess.PIPE, text=True,
        )
        for i in range(args.workers)
    ]
    counts = []
    for worker in workers:
        output, _ = worker.communicate()
        assert worker.returncode == 0
        counts.append(int(output))
    print(f"Python {platform.python_version()}, {len(workers)} workers, "
          f"{sum(counts)} grid tuples, {orbit_checks()} exact orbit checks, "
          f"{monotonic()-start:.3f}s elapsed")


if __name__ == "__main__":
    main()
