"""Rational lower bound for the noncentral contribution of a flat core."""
from fractions import Fraction as Q
from pathlib import Path
import json


def partial(x, degree):
    term = total = Q(1)
    for k in range(1, degree+1):
        term *= x/k
        total += term
    return total, term


def upper_exp(x):
    degree = 40
    total, term = partial(x, degree)
    first = term*x/(degree+1)
    ratio = x/(degree+2)
    assert 0 <= ratio < 1
    return total+first/(1-ratio)


def round_up(x, units):
    scaled = x*units
    return Q(-(-scaled.numerator//scaled.denominator), units)


def main():
    # J is the integral of p_{G(t)}(t) over (0,1), with G(t)=exp(-k(t)).
    # On each panel H(t)=t² exp(k(t))/2 has its maximum at an endpoint,
    # since its only positive critical point 1/8 is a minimum.
    count = 4096
    step = Q(39, 40*count)
    units = 10**12
    normalization = Q(3988, 10000)
    assert normalization**2 * 2*Q(22, 7) < 1
    ts = [Q(1, 40)+j*step for j in range(count+1)]
    ks = [(1/t-1)/4 for t in ts]
    upper_h = [round_up(t*t*upper_exp(k)/2, units) for t,k in zip(ts, ks)]
    sum_int = 0
    for j in range(count):
        exponent = ks[j+1]/2-max(upper_h[j], upper_h[j+1])
        lower = partial(exponent, 35)[0]
        assert lower > 0
        scaled = lower*units
        sum_int += scaled.numerator//scaled.denominator
    answer = normalization*step*Q(sum_int, units)
    assert answer > Q(59,100)
    output = {
        "panel_count": count, "start": "1/40", "integer_sum": sum_int,
        "exact_lower_bound": str(answer), "decimal_lower_bound": float(answer),
        "strict_lower_bound": "59/100",
    }
    Path("scaled-certificate.json").write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
