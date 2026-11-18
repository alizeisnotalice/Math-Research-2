from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from fractions import Fraction
from math import ceil, factorial, isqrt
from pathlib import Path
from time import monotonic, sleep
import argparse
import json
import os
import platform
import subprocess

import mpmath as mp


ROOT = Path(__file__).parent
RADICAL_BITS = 80
TERMS = 22


def utc():
    return datetime.now(timezone.utc).isoformat()


def sqrt_upper(value):
    scale = 1 << RADICAL_BITS
    scaled = (value.numerator << (2 * RADICAL_BITS)) // value.denominator
    return Fraction(isqrt(scaled) + 1, scale)


def bounds(k, eps, sigma_ratio):
    raw = Fraction(0)
    improved = Fraction(0)
    odds_reserved = sigma_ratio * sigma_ratio
    step = eps / 2
    factor = Fraction(1)
    falling = 1
    for ell in range(1, min(k // 2, TERMS) + 1):
        falling *= (k - 2 * ell + 2) * (k - 2 * ell + 1)
        factor *= step / (ell + 1)
        full = sqrt_upper(Fraction(falling))
        raw += factor * full
        if ell % 2:
            constant_weight = factorial(2 * ell) if ell > 1 else 0
            radicand = falling - (falling - constant_weight) * odds_reserved**ell
            improved += factor * sqrt_upper(radicand)
        else:
            improved += factor * full
    if k // 2 > TERMS:
        # All omitted positive summands are bounded by a geometric majorant.
        tail_ratio = step * k / (TERMS + 2)
        assert 0 < tail_ratio < 1
        tail = factor * full * tail_ratio / (1 - tail_ratio)
        raw += tail
        improved += tail
    mass_two = (eps * sigma_ratio)**2 / 8
    base = min(Fraction(1), raw * raw)
    refined = min(base, mass_two + improved * improved)
    return base, refined


def rational_profile(eps, r, limit):
    moments = [[Fraction(0), Fraction(0)] for _ in range(2)]
    single = [Fraction(0), Fraction(0)]
    argmax = [None, None]
    minimum_k = [None, None]
    for k in range(limit + 1):
        if k < 2:
            caps = [Fraction(0), Fraction(0)]
        elif k == 2:
            caps = [bounds(k, eps, r)[0], (eps*r)**2/8]
        else:
            caps = bounds(k, eps, r)
        for which, cap in enumerate(caps):
            tail = 1 - cap
            moments[which][0] += tail
            moments[which][1] += (2*k + 1) * tail
            if (k + 1) * tail > single[which]:
                single[which] = (k + 1) * tail
                argmax[which] = k
            if cap == 1 and minimum_k[which] is None:
                minimum_k[which] = k
        if minimum_k[1] is not None:
            break
    # Convert all exact rational lower sums downward on a fixed decimal grid.
    def downward(number):
        denominator = 10**16
        scaled = number.numerator * denominator // number.denominator
        return f"{scaled//denominator}.{scaled%denominator:016d}"

    return dict(epsilon=str(eps), sigma_ratio=str(r), evaluated_through=k,
                first_saturation=minimum_k, optimum_single_cutoff=argmax,
                single_lower=[downward(v) for v in single],
                energy_k_lower=[downward(v[0]) for v in moments],
                energy_k2_lower=[downward(v[1]) for v in moments])


def independent_profile(eps, r, limit):
    mp.mp.dps = 49
    eps = mp.mpf(eps.numerator) / eps.denominator
    r = mp.mpf(r.numerator) / r.denominator
    energy = [mp.mpf(0), mp.mpf(0)]
    for k in range(limit + 1):
        if k < 2:
            caps = [mp.mpf(0), mp.mpf(0)]
        elif k == 2:
            caps = [eps**2 / 8, (eps*r)**2 / 8]
        else:
            raw = mp.mpf(0)
            improved = mp.mpf(0)
            k_factorial = factorial(k)
            # Terms through 32 differ from the full sum by less than 1e-24
            # for this finite grid. The exact rational path handles the tail.
            for ell in range(1, min(k//2, 32) + 1):
                falling = k_factorial // factorial(k - 2*ell)
                factor = (eps/2)**ell / factorial(ell + 1)
                raw += factor * mp.sqrt(falling)
                if ell % 2:
                    constant = factorial(2*ell) if ell > 1 else 0
                    improved += factor * mp.sqrt(falling * (1-r**(2*ell)) + constant*r**(2*ell))
                else:
                    improved += factor * mp.sqrt(falling)
            caps = [min(1, raw**2), min(1, raw**2, (eps*r)**2/8 + improved**2)]
        for which, cap in enumerate(caps):
            energy[which] += 1-cap
        if caps[1] == 1:
            break
    return energy


def worker(payload):
    shard, duration, max_profiles = payload
    start_utc = utc()
    start = monotonic()
    result = ROOT / f'shard-{shard:02d}.json'
    rows = []
    for profile in range(max_profiles):
        index = 12*profile + shard
        shift = 6 + index % 4
        numerator = 8 + (index*19 // 4) % 17
        eps = Fraction(numerator, 8 * (1 << shift))
        r = Fraction(((index % 31) * 7) % 32, 32)
        cap_k = min(ceil(6 / eps), 1200)
        row = rational_profile(eps, r, cap_k)
        if profile < 1 or profile % 9 == 0:
            reference = independent_profile(eps, r, row['evaluated_through'])
            gaps = [reference[which] - mp.mpf(row['energy_k_lower'][which]) for which in range(2)]
            assert all(gap >= -mp.mpf('1e-22') and gap < mp.mpf('1e-14') for gap in gaps), gaps
            row['independent_mpmath_gap'] = [mp.nstr(gap, 10) for gap in gaps]
        rows.append(row)
        if monotonic() - start >= duration:
            break
    out = dict(pid=os.getpid(), started_utc=start_utc, ended_utc=utc(),
               duration_seconds=monotonic()-start, shard=shard, rows=rows)
    result.write_text(json.dumps(out, indent=2) + '\n')
    return out


def monitor(pids):
    ps = subprocess.run(['ps', '-p', ','.join(map(str, pids)),
                         '-o', 'pid=,%cpu=,state=,etime='], capture_output=True, text=True)
    top = subprocess.run(['top', '-l', '2', '-s', '1', '-n', '0'],
                         capture_output=True, text=True)
    load = [line for line in top.stdout.splitlines() if 'CPU usage:' in line]
    return dict(utc=utc(), workers=ps.stdout.splitlines(), cpu_usage=load[-1:])


def final_constants():
    mp.mp.dps = 60
    result = []
    for r in [0, mp.mpf('.5'), mp.mpf('.9'), 1]:
        def series(x):
            return sum(x**j * (mp.sqrt(1-r**(2*j)) if j % 2 else 1)
                       / mp.factorial(j+1) for j in range(1, 96))

        cutoff = mp.findroot(lambda x: series(x) - 1, (1, 3))
        best = mp.findroot(lambda x: 1-series(x)**2 - 2*x*series(x)*mp.diff(series, x),
                           (.5, cutoff*.95))
        scale1 = 2 * mp.quad(lambda x: 1-series(x)**2, [0, cutoff])
        scale2 = 8 * mp.quad(lambda x: x*(1-series(x)**2), [0, cutoff])
        result.append(dict(r=str(r), x_cutoff=mp.nstr(cutoff, 42),
                           x_best_single=mp.nstr(best, 42),
                           epsilon_times_single_lower=mp.nstr(2*best*(1-series(best)**2), 42),
                           epsilon_times_energy_k_lower=mp.nstr(scale1, 42),
                           epsilon2_times_energy_k2_lower=mp.nstr(scale2, 42)))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--duration', type=float, default=55)
    parser.add_argument('--profiles', type=int, default=240)
    args = parser.parse_args()
    start = utc()
    with ProcessPoolExecutor(max_workers=12) as executor:
        futures = [executor.submit(worker, (i, args.duration, args.profiles)) for i in range(12)]
        processes = executor._processes
        pids = list(processes)
        with (ROOT / "telemetry.jsonl").open("w") as file:
            print(json.dumps(dict(utc=utc(), parent=os.getpid(), pids=pids,
                                  platform=platform.platform(), python=platform.python_version(),
                                  mpmath=mp.__version__, cpu_count=os.cpu_count())), file=file, flush=True)
            while not all(j.done() for j in futures):
                print(json.dumps(monitor(pids)), file=file, flush=True)
                sleep(1.1)
        rows = [j.result() for j in futures]
    after = subprocess.run(['ps','-p',','.join(map(str,pids)), '-o','pid='],
                           capture_output=True, text=True)
    result = dict(started_utc=start, workers=[{k:v for k,v in row.items() if k!='rows'}
                                            for row in rows],
                  end_workers_utc=utc(), rows=sum(len(row['rows']) for row in rows),
                  independent_checks=sum('independent_mpmath_gap' in element
                                         for row in rows for element in row['rows']),
                  worker_exit_codes={pid:proc.exitcode for pid,proc in processes.items()},
                  processes_after=after.stdout.splitlines(), scaled_constants=final_constants())
    (ROOT/'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
