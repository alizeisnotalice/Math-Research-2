"""Independent adaptive checks for the largest quartic midpoint-scan ratios."""
import argparse
import json
import math
from pathlib import Path
import platform
import sys
import warnings

import mpmath as mp
import scipy
from scipy.integrate import IntegrationWarning, quad


def parts(hit):
    radius = hit["radius"]
    split = max(-radius, min(radius, hit["shift"]))
    return [(-radius, split, hit["coefficient"], 4), (split, radius, hit["coefficient"], 4)]


def adaptive_float64(hit):
    x, shift, epsilon, variance = (hit[k] for k in ("center", "shift", "epsilon", "variance"))
    total = 0.0
    error = 0.0
    for left, right, scale, power in parts(hit):
        if left == right:
            continue

        def integrand(u):
            v = scale * abs(u) ** power + epsilon
            log_ratio = 0.5 * (math.log((variance + epsilon) / (2 * math.pi * v))
                               + (x + shift) ** 2 / (variance + epsilon)
                               - (x + u) ** 2 - (shift - u) ** 2 / v)
            return math.exp(log_ratio)

        value, err = quad(integrand, left, right, epsabs=1e-12, epsrel=1e-12, limit=300)
        total += value
        error += err
    return total, error


def adaptive_mp(hit, dps):
    with mp.workdps(dps):
        convert = lambda key: mp.mpf(str(hit[key]))
        x, shift, epsilon, variance = map(convert, ("center", "shift", "epsilon", "variance"))
        total = mp.mpf(0)
        for left, right, scale, power in parts(hit):
            if left == right:
                continue
            left, right, scale = (mp.mpf(str(z)) for z in (left, right, scale))

            def integrand(u):
                v = scale * abs(u) ** power + epsilon
                log_ratio = ((x + shift) ** 2 / (variance + epsilon) - (x + u) ** 2 - (shift - u) ** 2 / v) / 2
                return mp.sqrt((variance + epsilon) / (2 * mp.pi * v)) * mp.exp(log_ratio)

            total += mp.quad(integrand, [left, (left + right) / 2, right])
        return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("json_input", type=Path)
    parser.add_argument("json_output", type=Path)
    args = parser.parse_args()
    grid = json.loads(args.json_input.read_text())
    candidates = grid["bestByRegime"]
    results = {}
    warnings.simplefilter("error", IntegrationWarning)
    for regime, hit in candidates.items():
        value, error = adaptive_float64(hit)
        high_precision = adaptive_mp(hit, 70)
        lower_precision = adaptive_mp(hit, 50)
        with mp.workdps(70):
            mp_error = abs(high_precision - lower_precision)
        if abs(value - float(high_precision)) > 1e-10 * max(1.0, abs(value)):
            raise RuntimeError(f"Adaptive quadratures disagree for {regime}")
        results[regime] = {"parameters": hit, "scipyAdaptive": value,
                           "scipyEstimatedError": error,
                           "mpmath70Digits": mp.nstr(high_precision, 65),
                           "mpmath50Vs70AbsDifference": mp.nstr(mp_error, 12),
                           "gpuVsCPUAbsDifference": abs(value - hit["ratio"])}
        print(f"{regime}: metal={hit['ratio']:.8g} scipy={value:.14g}, disagreement={abs(value-hit['ratio']):.3g}")
    report = {"python": sys.version, "scipy": scipy.__version__, "mpmath": mp.__version__,
              "architecture": platform.machine(), "comparison": results}
    args.json_output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
