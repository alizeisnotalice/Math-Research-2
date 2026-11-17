#!/usr/bin/env python3
"""Independent double-precision audit of the R018 local saddle formula.

Unlike round018_log_contact.py, this integrates directly in w=t/t_* with
SciPy's adaptive Gauss--Kronrod routine.  It evaluates exp(g(t_*)) times the
local t-integral, avoiding the exponentially small unscaled integral.
"""
from __future__ import annotations

import math
import platform

import scipy
from scipy.integrate import quad

X0 = 3.0
C = 1.0
T0 = 4.0


def density(x: float) -> float:
    return math.exp(-x * x / 2.0) / math.sqrt(2.0 * math.pi)


def scaled_local_integral(beta: float, xi: float) -> tuple[float, float, float, float]:
    """Return exp(g*) I, its quadrature error, leading term, and ratio.

    I is exactly the local integral over 0 < |u| < exp(-T0), so that
    t=log(1/|u|) ranges over [T0,infinity).  In w=t/t*,
    exp(g*) I has integrand t* h(t*w) exp(-t* Phi(w)).
    """
    tstar = (C * beta * xi * xi / 2.0) ** (1.0 / (beta + 1.0))
    gstar = (beta + 1.0) * tstar / beta
    lower = T0 / tstar

    def fn(w: float) -> float:
        t = tstar * w
        u = math.exp(-t)
        phi = w + w ** (-beta) / beta - 1.0 - 1.0 / beta
        return (
            tstar
            * (density(X0 + u) + density(X0 - u))
            * math.cos(xi * u)
            * math.exp(-tstar * phi)
        )

    cuts = [lower]
    for a in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
        if a > lower:
            cuts.append(a)
    cuts.append(math.inf)
    value = 0.0
    error = 0.0
    for left, right in zip(cuts, cuts[1:]):
        piece, piece_error = quad(fn, left, right, epsabs=2e-13, epsrel=2e-13, limit=400)
        value += piece
        error += piece_error

    leading = 2.0 * density(X0) * math.sqrt(2.0 * math.pi * tstar / (beta + 1.0))
    return value, error, leading, value / leading


def main() -> None:
    print(f"python={platform.python_version()} scipy={scipy.__version__}")
    print("method=direct w=t/t_star; scipy.integrate.quad adaptive Gauss-Kronrod")
    print("beta xi t_star scaled_integral error_bound leading ratio")
    for beta, xi in ((1.0, 100.0), (2.0, 500.0), (0.5, 100.0)):
        tstar = (C * beta * xi * xi / 2.0) ** (1.0 / (beta + 1.0))
        value, error, leading, ratio = scaled_local_integral(beta, xi)
        print(
            f"{beta:g} {xi:g} {tstar:.16g} {value:.16g} {error:.2e} "
            f"{leading:.16g} {ratio:.16g}"
        )


if __name__ == "__main__":
    main()
