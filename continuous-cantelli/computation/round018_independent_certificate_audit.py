"""Independent, non-mutating audit of round 018's finite-frequency certificate.

This file deliberately does not import or execute round018_log_contact.py, because
that program writes its result JSON.  It rechecks the exact comparisons and the
displayed saddle diagnostic formulas against the checked-in JSON.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib
import json

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "computation" / "round018_log_contact.py"
RESULT = ROOT / "computation" / "round018-log-results.json"


def exact_audit() -> None:
    # e > 2 follows from its positive series tail; e < 3 follows by bounding
    # sum_{n>=3} 1/n! < sum_{n>=3} 1/(2*3^(n-2)) = 1/4.
    assert Q(5, 2) + Q(1, 4) < 3
    assert 2**100 > 200
    # e^5 > 100 and cos(1)>1/2 by, respectively, a positive exponential
    # partial sum and an alternating cosine lower partial sum.
    assert sum(Q(5**j, factorial(j)) for j in range(7)) > 100
    assert Q(1) - Q(1, 2) + Q(1, 24) - Q(1, 720) > Q(1, 2)
    # pi < 22/7 gives sqrt(2*pi) < 3, since 44/7 < 9.
    assert Q(44, 7) < 9
    # On the two annuli total length is 2e^-100; the 1/2 cosine lower bound
    # cancels the factor 2.  Thus the analytic lower bound is
    # 3^-9 e^-151 > 3^-160 (using e < 3).
    positive = Q(1, 3**160)
    negative = Q(1, 2**1000)
    target = Q(1, 2**5000)
    assert Q(10000, 198) < 51
    assert positive - negative > positive / 2 > target
    # e^4 < 3^4 < 100 puts |u|<=1/100 in the stated local model;
    # log(100)<5 follows from e^5>100.
    assert 3**4 < 100


def scaled_local(beta: str, frequency: int, dps: int = 100) -> mp.mpf:
    """Compute exp(g(t_*)) J(xi), independently matching formula (14)."""
    mp.mp.dps = dps
    b = mp.mpf(beta)
    xi = mp.mpf(frequency)
    tstar = (b * xi**2 / 2) ** (1 / (b + 1))
    root = mp.sqrt(tstar)
    lower = (mp.mpf(4) - tstar) / root
    normal = lambda x: mp.exp(-x * x / 2) / mp.sqrt(2 * mp.pi)

    def integrand(y: mp.mpf) -> mp.mpf:
        t = tstar + root * y
        w = t / tstar
        g_minus_gstar = tstar * (w + w ** (-b) / b - 1 - 1 / b)
        u = mp.exp(-t)
        return root * (normal(3 + u) + normal(3 - u)) * mp.cos(xi * u) * mp.exp(-g_minus_gstar)

    knots = [lower]
    knots += [mp.mpf(k) for k in (-40, -20, -12, -8, -4, 0, 4, 8, 12, 20, 40) if k > lower]
    knots.append(mp.inf)
    return mp.quad(integrand, knots)


def numerical_audit() -> None:
    record = json.loads(RESULT.read_text())
    by_key = {(row["beta"], row["frequency"], row["dps"]): row for row in record["diagnostics"]}
    # Every stored 80-digit integral is recomputed independently at 100 dps.
    for beta in ("0.5", "1", "2"):
        for xi in (20, 100, 500):
            got = scaled_local(beta, xi)
            stored = mp.mpf(by_key[(beta, xi, 80)]["scaled_local_characteristic"])
            assert abs(got - stored) < mp.mpf("1e-70"), (beta, xi, got, stored)


if __name__ == "__main__":
    exact_audit()
    numerical_audit()
    result = json.loads(RESULT.read_text())
    assert result["script_sha256"] == hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    print("round018 independent audit passed")
    print("source_sha256", result["script_sha256"])
    print("mpmath", mp.__version__)
