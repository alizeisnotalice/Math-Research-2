"""Independent mpmath quadrature and exploratory fit of log absolute ratio error."""

import mpmath as mp
from approximants import aptekarev
from sondow import case


def integral(n):
    def kernel(x, y):
        if x <= 0 or y <= 0 or x >= 1 or y >= 1:
            return 0
        u = x * y
        return (x * (1 - x) * y * (1 - y)) ** n / ((1 - u) * (-mp.log(u)))

    return mp.quad(lambda x: mp.quad(lambda y: kernel(x, y), [0, 1]), [0, 1])


for precision in (22, 32):
    with mp.workdps(precision):
        for n in (1, 2, 3):
            row = case(n)
            value = integral(n)
            logs = mp.fsum(v * mp.log(p) for p, v in row["prime_exponents"].items())
            identity = (row["D"] * mp.euler + logs - row["dA"]) / row["d"]
            residual = abs(value - identity)
            if residual > mp.mpf(10) ** (7 - precision):
                raise ArithmeticError((n, residual))
            print(
                precision,
                n,
                "I",
                mp.nstr(value, 18),
                "absolute residual",
                mp.nstr(residual, 4),
            )
with mp.workdps(80):
    for first, last in [(16, 64), (64, 144)]:
        ns = list(range(first, last + 1, 8))
        xs = mp.matrix([[1, mp.sqrt(n)] for n in ns])
        pairs = [aptekarev(n) for n in ns]
        ys = mp.matrix([mp.log(abs(mp.euler - mp.mpf(p) / q)) for p, q in pairs])
        beta, res = mp.qr_solve(xs, ys)
        print(
            "fit",
            first,
            last,
            "slope",
            mp.nstr(beta[1], 12),
            "target",
            mp.nstr(-2 * mp.sqrt(2), 12),
        )
