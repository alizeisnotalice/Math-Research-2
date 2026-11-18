#!/usr/bin/env python3
"""SciPy adaptive double integration of independently evaluated trapezoid segments."""
import argparse
import json
import math
import re
import statistics
import struct
import sys
from pathlib import Path

import scipy
from scipy.integrate import quad
from scipy.special import erf, logsumexp

LOG_NORM = -0.5 * math.log(2 * math.pi)
LOG_WINDOW_MASS = math.log(erf(1 / math.sqrt(2)))


def float32(x):
    return struct.unpack("=f", struct.pack("=f", x))[0]


def integral(point, unit):
    rho, power, a, delta, gamma, lam = (
        float32(point[k]) for k in ("rho", "p", "a", "delta", "gamma", "lambda")
    )
    q = lam**2 / 2
    heat = q * (delta + gamma)
    assert heat > 0
    logs = [-q * delta + (LOG_WINDOW_MASS if unit else 0.0)]
    for j in range(12):
        center = 4 * rho**j
        diff = center - lam
        log_width = -math.log(8) - a * center**power
        width = math.exp(log_width)
        left, right = -1., 1.
        if unit and width:
            left = max(left, (-1 - diff) / width)
            right = min(right, (1 - diff) / width)
        elif unit:
            if abs(diff) > 1:
                continue
            if diff == 1:
                right = 0.
            if diff == -1:
                left = 0.
        if right <= left:
            continue

        for start, end, intercept, beta in ((-1., -.5, 2., 2.),
                                             (-.5, .5, 1., 0.),
                                             (.5, 1., 2., -2.)):
            lo, hi = max(left, start), min(right, end)
            if hi <= lo:
                continue
            # Orient steep slopes into the integration interval from its inside edge.
            x0, x1 = (hi, lo) if beta > 0 else (lo, hi)
            dx = x1 - x0
            peak = intercept + beta * x0
            fall = -beta * dx
            d0 = diff + width * x0
            span = width * dx
            slope = heat * fall + d0 * span
            curvature = span * span
            vertex = max(0., min(1., -slope / curvature)) if curvature else float(slope < 0)
            offset = -slope * vertex - 0.5 * curvature * vertex**2
            derivative = -slope - curvature * vertex
            value = 0.

            for direction, length in ((-1., vertex), (1., 1. - vertex)):
                if length <= 0:
                    continue
                scale = max(1., -direction * derivative, math.sqrt(curvature))
                extent = min(scale * length, 100.)

                def scaled(z):
                    displacement = direction * z / scale
                    t = vertex + displacement
                    v = heat * max(0., peak - fall * t)
                    exponent = derivative * displacement - 0.5 * curvature * displacement**2
                    return math.exp(exponent) * (-math.expm1(-v))

                part, error = quad(scaled, 0, extent, epsabs=5e-14, epsrel=5e-13, limit=120)
                assert part >= 0 and error < 1e-9 * max(1., part)
                value += part / scale
            if value:
                logs.append(log_width + q * gamma + LOG_NORM - 0.5*d0*d0
                            - heat * (1 - peak) + offset + math.log(abs(dx) * value))
    return float(logsumexp(logs))


def audit(data):
    checked = {}
    for key, row in data["summaries"].items():
        hit = row["best"]
        if hit is None:
            continue
        point = hit["p"]
        lam2 = float32(point["lambda"])**2
        exact = [integral(point, unit) for unit in (False, True)]
        gpu = [float32(hit[k]) for k in ("logFull", "logUnit")]
        differences = [a-b for a,b in zip(gpu, exact)]
        checked[key] = dict(point=point, index=hit["index"], salt=hit["salt"],
                            cpu_full=exact[0], cpu_unit=exact[1], gpu_full=gpu[0],
                            gpu_unit=gpu[1], gpu_minus_cpu=differences,
                            normalized_difference=[d/lam2 for d in differences],
                            signs_agree=[math.copysign(1,a)==math.copysign(1,b) and (a==0)==(b==0)
                                         for a,b in zip(gpu,exact)])
    return checked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gpu_json", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--ioreg", type=Path)
    args = parser.parse_args()
    data = json.loads(args.gpu_json.read_text())
    checked = audit(data)
    controls = {}
    for lam in (256., 320., 512.):
        point = dict(rho=4., p=2., a=1/32, delta=1/4, gamma=1/8, **{"lambda":lam})
        controls[str(lam)] = [integral(point, unit) for unit in (False,True)]
    result = dict(source=args.gpu_json.name, python=sys.version, scipy=scipy.__version__,
                  cpu_check=checked, fixed_profile_log_full_and_unit=controls,
                  max_normalized_difference=max(
                      abs(d) for row in checked.values() for d in row["normalized_difference"]))
    if args.ioreg is not None:
        used = [int(m.group(1)) for line in args.ioreg.read_text().splitlines()
                if line.startswith("during ") and (m := re.search(r'"Device Utilization %"=(\d+)', line))]
        result["ioreg_during"] = dict(count=len(used), minimum=min(used),
                                       median=statistics.median(used), maximum=max(used))
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("max normalized abs discrepancy", result["max_normalized_difference"])
    for key,row in checked.items():
        print(key, row["cpu_full"], row["cpu_unit"], row["signs_agree"])
    print("control profiles", controls)


if __name__ == "__main__":
    main()
