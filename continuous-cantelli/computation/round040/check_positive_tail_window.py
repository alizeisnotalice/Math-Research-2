"""Check the mark identity and calibrate two local error mechanisms."""

import json
import sys

import mpmath as mp
import numpy as np
import scipy
from scipy.special import roots_hermitenorm, roots_legendre
import sympy as sp


def exact_checks():
    t, s, v = sp.symbols("t s v", positive=True)
    z = t + sp.I * s
    assert sp.expand(-z**2 * v / (2*t**2) + v/2 + sp.I*s*v/t - s**2*v/(2*t**2)) == 0
    y, a, b = sp.symbols("y a b", real=True)
    primitive = sp.Function("P")
    derivative = sp.diff(primitive(y-a) - primitive(y-b), y)
    assert derivative == sp.Subs(sp.Derivative(primitive(y), y), y, y-a) - sp.Subs(
        sp.Derivative(primitive(y), y), y, y-b
    )


def kernel(y):
    return (
        np.sinc(y / np.pi)**4
        + np.sinc((y + np.pi / 2) / np.pi)**4
        + np.sinc((y - np.pi / 2) / np.pi)**4
    )


def mark_comparison():
    legendre_x, legendre_w = roots_legendre(72)
    normal_x, normal_w = roots_hermitenorm(48)
    normal_w /= np.sqrt(2 * np.pi)
    y = np.array([0.0, 0.1, 0.7, 1.6, 3.1, 6.3, 9.2, 15.7, 32.95, 64.0])
    cutoffs = [0.1, 1.0, 8.0, 40.0]
    results = []
    for L in (8, 16, 32, 64):
        worst = 0.0
        for cutoff in cutoffs:
            mark = (legendre_x + 1) * cutoff / 2
            mark_weight = legendre_w * cutoff / 4 * np.exp(-mark / 2)
            mark_weight /= -np.expm1(-cutoff / 2)
            np.testing.assert_allclose(mark_weight.sum(), 1.0, rtol=1e-12)
            w = (-mark[:, None] + 1j * np.sqrt(mark[:, None]) * normal_x[None, :]) / L
            shifted = kernel(y[:, None, None] + w[None, :, :])
            values = np.einsum("v,g,yvg->y", mark_weight, normal_w, shifted)
            largest = np.max(np.abs(values / kernel(y) - 1))
            worst = max(worst, float(largest))
        results.append({"L": L, "max_relative_error_on_fixed_grid": worst, "L_times_error": L*worst})
    return results


def compact_calibration():
    A = mp.exp(mp.mpf(1)/2)
    rows = []
    for t in (16, 64, 256, 1024):
        s = mp.mpf(1)
        z = t + 1j*s
        compact_mass = mp.exp(t**2 / (2*(1+t**2))) / mp.sqrt(1+t**2)
        complex_compact = mp.exp(z**2 / (2*(1+z**2))) / mp.sqrt(1+z**2)
        diff = t**2 / z**2 * complex_compact - compact_mass
        rows.append({"t": t, "t_squared_times_imag_residual_without_gaussian": mp.nstr(t**2 * diff.imag, 20)})
    return rows, mp.nstr(-3*A, 20)


def critical_window():
    # The positive profile B(x)=A/x*(1+sin(x*x)/2) with A=1.
    d = mp.pi/2
    rows = []
    for n in (10, 100, 10000):
        t = mp.sqrt(2*mp.pi*n)
        ratio = mp.quad(
            lambda u: mp.sin(2*u + u*u/t**2) / (2*(1+u/t**2)),
            [0, d],
        )
        rows.append({"n": n, "t_squared_times_window_residual": mp.nstr(ratio, 20)})
    return rows


def main():
    exact_checks()
    mark = mark_comparison()
    with mp.workdps(80):
        compact, compact_limit = compact_calibration()
        result = {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
            "mark_quadrature": mark,
            "compact_model": compact,
            "compact_model_limit": compact_limit,
            "synthetic_critical_window": critical_window(),
            "synthetic_critical_window_limit": "0.5",
        }
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
