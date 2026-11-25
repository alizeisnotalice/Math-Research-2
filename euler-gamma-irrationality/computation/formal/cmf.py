"""Derive row-state contiguous shifts at z=1; compare two dual trajectories."""

import json
import sys
from pathlib import Path
from math import gcd
import sympy as s

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from intervals import certify


def shift(a, b, c, sign):
    # theta^3 expressed in (1,theta,theta^2), from the Meijer differential equation.
    x, y = a
    i, j, k = b
    theta = s.Matrix(
        [
            [0, 0, i * j * k - (1 - x) * (1 - y)],
            [1, 0, -i * j - i * k - j * k - 2 + x + y],
            [0, 1, i + j + k - 1],
        ]
    )
    return c * s.eye(3) + sign * theta


def trajectory(j, extra):
    a, b = [-extra * j, 0], [j] * 3
    M = s.eye(3)
    if extra:
        M = M * shift(a, b, 1 - a[0], 1)
        a[0] -= 1
    for k in range(3):
        M = M * shift(a, b, b[k], -1)
        b[k] += 1
    return M.applyfunc(s.factor)


def run(extra):
    j = s.symbols("j", integer=True, nonnegative=True)
    M = trajectory(j, extra)
    determinant = s.factor(M.det())
    D = M.inv().T.applyfunc(s.factor)
    scale = s.factor(s.lcm([s.denom(z) for z in D]))
    D = (scale * D).applyfunc(s.factor)
    assert (M * D.T - scale * s.eye(3)).applyfunc(s.factor) == s.zeros(3)
    u, v = s.Matrix([[1, 1, 0]]), s.Matrix([[1, 1, 1]])
    out = []
    for n in range(1, 33):
        step = D.subs(j, n - 1)
        u, v = u * step, v * step
        P, Q = int(u[0]), int(v[0])
        assert Q
        g = gcd(P, Q)
        p, q = P // g, Q // g
        if q < 0:
            p, q = -p, -q
        bits = q.bit_length() + 80
        out.append(
            dict(
                n=n,
                P=P,
                Q=Q,
                gcd=g,
                p=p,
                q=q,
                form=certify(q, -p, {}, start=bits, limit=max(bits, 512)),
            )
        )
    return dict(
        extra=extra, det=str(determinant), scale=str(scale), D=str(D), pairs=out
    )


if __name__ == "__main__":
    print(json.dumps([run(0), run(1)], indent=2))
