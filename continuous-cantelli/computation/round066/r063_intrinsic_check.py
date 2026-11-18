"""Deterministic quadratures of the two intrinsic-scale multipliers."""
import json
import platform
from pathlib import Path

import mpmath as m

m.mp.dps = 50
c = m.mpf(1)
gaussian_normalization = 1 / m.sqrt(2 * m.pi)

def density(x, variance=1):
    return gaussian_normalization * m.exp(-x**2 / (2 * variance)) / m.sqrt(variance)

def density_d4(x, variance):
    z2 = x**2 / variance
    return (z2**2 - 6*z2 + 3) * density(x, variance) / variance**2

def num(x):
    return m.nstr(x, 24)

amplitude = m.mpf(2) / 5
def harmonic(t):
    return m.cos(t) * m.exp(-(c + amplitude * m.cos(t)) / 2)
first_harmonic = m.quad(harmonic, [0, m.pi, 2 * m.pi]) / (2 * m.pi)
exact_harmonic = -m.exp(-c / 2) * m.besseli(1, amplitude / 2)
predicted = -gaussian_normalization * m.exp(-(1+c) / 2) / 2
rows = []
for n in (10, 30, 100, 300):
    r = 1 / (2 * m.pi * n)
    def response(t):
        factor = 1 + r*t
        u = r * factor
        baseline = 1 + c * factor**4
        variance = baseline + factor**4 * u**2 * m.cos(1/u)
        return density(u) * (density_d4(-t, variance) - density_d4(-t, baseline)) / r**2
    approx = m.quad(response, list(range(-20, 21, 2)))
    rows.append({"n": n, "r": num(r), "normalized_response": num(approx),
                 "absolute_error": num(abs(approx - predicted))})
assert abs(first_harmonic - exact_harmonic) < m.mpf("1e-45")
assert m.mpf(rows[-1]["absolute_error"]) < m.mpf("2e-5")
print(json.dumps({"python": platform.python_version(), "mpmath": m.__version__,
                  "decimal_precision": m.mp.dps, "randomness": "none",
                  "fixed_variance_harmonic": num(first_harmonic),
                  "bessel_absolute_error": num(abs(first_harmonic-exact_harmonic)),
                  "q": 4, "c": 1, "k": 4, "omega": 1, "phase": 0,
                  "H": "cos", "eta": "u^2", "g": "standard normal",
                  "s": 0, "integration_interval": [-20, 20],
                  "segments": 20, "limiting_normalized_response": num(predicted),
                  "rows": rows}, indent=2))
