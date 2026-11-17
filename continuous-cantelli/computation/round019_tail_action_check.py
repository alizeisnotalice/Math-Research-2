#!/usr/bin/env python3
"""Symbolic and numerical audit for the R019 critical two-tail estimate."""
from __future__ import annotations

import platform

import mpmath as mp
import sympy as sp


def constants(beta: mp.mpf, d: mp.mpf, eps: mp.mpf, C: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    p = 2 / (beta + 1)
    aeps = (beta + 1) / (2 * beta) * (1 - eps) ** (beta / (beta + 1)) * (beta * d) ** (1 / (beta + 1))
    kloc = (beta + 1) / beta * (beta * C / 2) ** (1 / (beta + 1))
    return p, aeps, kloc


def equality_tail(beta: mp.mpf, d: mp.mpf, q: mp.mpf, radius: mp.mpf) -> mp.mpf:
    """Integral over |x|>=radius for v(x)=d|x|^(-2 beta)."""
    xstar = (beta * d * q * q) ** (1 / (2 * beta + 2))

    def side(x: mp.mpf) -> mp.mpf:
        return mp.exp(-x * x / 2 - d * q * q / (2 * x ** (2 * beta))) / mp.sqrt(2 * mp.pi)

    points = [radius]
    for point in (xstar / 2, xstar, 2 * xstar):
        if point > radius:
            points.append(point)
    points.append(mp.inf)
    return 2 * mp.quad(side, points)


def symbolic_check() -> None:
    b, d, q, e, s = sp.symbols("b d q e s", positive=True)
    phi = (1 - e) * s / 2 + d * q**2 / (2 * s**b)
    derivative = sp.diff(phi, s)
    stationary_relation = sp.solve(sp.Eq(derivative, 0), q**2)[0]
    minimum_from_relation = sp.simplify(phi.subs(q**2, stationary_relation))
    claimed_from_relation = sp.simplify(
        (b + 1) / (2 * b) * (1 - e) * s
    )
    assert sp.simplify(minimum_from_relation - claimed_from_relation) == 0
    print("sympy stationary q^2:", stationary_relation)
    print("sympy minimum difference after stationarity:", sp.simplify(minimum_from_relation - claimed_from_relation))


def main() -> None:
    mp.mp.dps = 70
    print(f"python={platform.python_version()} mpmath={mp.__version__} sympy={sp.__version__}")
    symbolic_check()
    print("beta d eps q R  tail_integral  full_line_bound  ratio  Aeps  Kloc")
    for beta, d, eps, q, radius in (
        (mp.mpf("1"), mp.mpf("4"), mp.mpf("0.1"), mp.mpf("20"), mp.mpf("1")),
        (mp.mpf("1.5"), mp.mpf("4"), mp.mpf("0.2"), mp.mpf("20"), mp.mpf("1")),
    ):
        p, aeps, kloc = constants(beta, d, eps, mp.mpf("1"))
        tail = equality_tail(beta, d, q, radius)
        bound = eps ** (-mp.mpf("0.5")) * mp.exp(-aeps * q**p)
        assert 0 < tail < bound
        print(
            mp.nstr(beta, 5), mp.nstr(d, 5), mp.nstr(eps, 5), mp.nstr(q, 5), mp.nstr(radius, 5),
            mp.nstr(tail, 18), mp.nstr(bound, 18), mp.nstr(tail / bound, 12),
            mp.nstr(aeps, 14), mp.nstr(kloc, 14),
        )

    # This demonstrates why d>2^beta C only permits a sufficiently small eps.
    beta = mp.mpf("1")
    d = mp.mpf("3")
    C = mp.mpf("1")
    for eps in (mp.mpf("0.1"), mp.mpf("0.5")):
        _, aeps, kloc = constants(beta, d, eps, C)
        print("fixed-eps comparison beta=1,d=3,C=1,eps=", eps, "Aeps/Kloc=", mp.nstr(aeps / kloc, 16))


if __name__ == "__main__":
    main()
