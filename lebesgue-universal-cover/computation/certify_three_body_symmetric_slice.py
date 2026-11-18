"""Arb certificate for the exact global minimum on the zero-twist slice."""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import platform

import flint
from flint import arb, ctx

ctx.prec = 160
pi = arb.pi()
r = arb(1) / 2
R3 = 1 / arb(3).sqrt()
R5 = 1 / (2 * (pi / 10).cos())
g = R5 * (pi / 5).sin()
c36 = R5 * (pi / 5).cos()
c72 = R5 * (2 * pi / 5).cos()


def q(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def interval(first, last):
    lo, hi = q(first), q(last)
    assert lo < hi
    return arb((lo + hi) / 2, (hi - lo) / 2)


def normal(angle):
    return angle.cos(), angle.sin()


def dot(x, y):
    return x[0] * y[0] + x[1] * y[1]


def det(x, y):
    return x[0] * y[1] - x[1] * y[0]


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def sub(x, y):
    return x[0] - y[0], x[1] - y[1]


def scale(s, x):
    return s * x[0], s * x[1]


def check_topology():
    # Strict inequalities on an open neighborhood of the supporting triple.
    a = interval('-.01216', '-.01215')
    b = interval('.02299', '.02301')
    right_B3 = R3 + a
    left_B3 = a - R3 / 2
    right_B5 = R5 + b
    right_upper_B5 = c72 + b
    w = (b - c36, -g)
    e18, e30, e54 = (normal(pi * k / 180) for k in (18, 30, 54))
    conditions = {
        'B3_right_vertex_beyond_B5_right': right_B3 - right_B5,
        'B5_right_vertex_beats_disk_through_18': right_B5 * e18[0] - r,
        'B3_vertex_loses_to_disk_by_30': r - right_B3 * e30[0],
        'B5_arc_loses_to_disk_by_30': -r - dot(w, e30),
        'B3_vertex_loses_to_B5_arc_by_30': 1 - dot((right_B3-w[0], g), e30),
        'B5_arc_support_decreases_at_30': w[0] * e30[1] + g * e30[0],
        'relative_B3_vertex_support_decreases_at_18': (right_B3-w[0]) * e18[1] - g * e18[0],
        'B5_arc_center_negative_x': -w[0],
        'B5_arc_center_negative_y': g,
        'relative_B3_arc_x_positive': left_B3 - w[0],
        'relative_B3_arc_y_negative': r - g,
        'B3_arc_loses_to_B5_arc_at_30': -dot((left_B3 - w[0], -r + g), e30),
        'B5_arc_beats_disk_at_54': r + dot(w, e54),
        'B5_vertex_beats_disk_at_54': dot((right_upper_B5, r), e54) - r,
    }
    for name, value in conditions.items():
        assert value > 0, (name, value)
    return {name: str(value) for name, value in conditions.items()}


def triple(a):
    c = R3 + a
    x = r / c
    s = (1 - x**2).sqrt()
    w = ((-r + g*s)/x, -g)
    b = w[0] + c36
    # Both gamma and delta have supporting equation w dot n = -1/2.
    n_gamma = (x, s)
    n_delta = sub(scale(-1/dot(w, w), w), n_gamma)
    return c, b, w, n_gamma, n_delta


def balanced_gradient(a):
    _, _, _, n_gamma, n_delta = triple(a)
    s = n_gamma[1]
    M = 1 + g - 3*n_delta[1]/2
    minus = (3*s/2 + g - 1, M-s/2)
    plus = (5*s/2 - g - 1, M+s/2)
    assert minus[1] < 0 < plus[1]
    weight_plus = -minus[1] / s
    assert 0 < weight_plus < 1
    gradient = minus[0] + weight_plus * (plus[0]-minus[0])
    return gradient, minus, plus, weight_plus


def triple_area(a):
    c, b, w, n_gamma, n_delta = triple(a)
    assert q('.02299') < b < q('.02301')
    assert 0 < n_delta[0] < n_gamma[0] < 1
    assert 0 < n_gamma[1] < n_delta[1] < 1
    gamma = (n_gamma[1]/n_gamma[0]).atan()
    delta = (n_delta[1]/n_delta[0]).atan()
    assert pi/10 < gamma < pi/6 < delta < 3*pi/10
    e54 = normal(3*pi/10)
    e150 = normal(5*pi/6)
    opposite_gamma = (-n_gamma[0], n_gamma[1])
    p0 = (c, arb(0))
    p5_up = (b+c72, r)
    p3_up = (a-R3/2, r)
    values = [
        det(p0, scale(r, n_gamma)),
        r**2 * (delta-gamma),
        det(scale(r, n_delta), add(w, n_delta)),
        det(w, sub(e54, n_delta)) + 3*pi/10-delta,
        det(p5_up, p3_up),
        det(p0, sub(opposite_gamma, e150)) + pi/6-gamma,
        det(add(p0, opposite_gamma), scale(r, opposite_gamma)),
        r**2 * gamma,
    ]
    return sum(values), b, gamma, delta


def main():
    topology = check_topology()
    a0 = '-.012154558715405409'
    a1 = '-.012154558715405408'
    zero_left, *_ = balanced_gradient(q(a0))
    zero_right, *_ = balanced_gradient(q(a1))
    assert zero_left < 0 < zero_right
    root_bracket = interval(a0, a1)
    _, lower_gradient, upper_gradient, weight_plus = balanced_gradient(root_bracket)
    area, b, gamma, delta = triple_area(root_bracket)
    area_lower, area_upper = '.8347809457752035', '.8347809457752037'
    assert q(area_lower) < area < q(area_upper), area
    assert q(area_upper) < q('.834780947305')
    source = Path(__file__)
    print(json.dumps({
        'python': platform.python_version(),
        'python_flint': flint.__version__,
        'precision_bits': ctx.prec,
        'script_sha256': sha256(source.read_bytes()).hexdigest(),
        'topology_sign_intervals': topology,
        'root_bracket': [a0, a1],
        'zero_left': str(zero_left),
        'zero_right': str(zero_right),
        'minus_gradient': list(map(str, lower_gradient)),
        'plus_gradient': list(map(str, upper_gradient)),
        'balanced_weight_plus': str(weight_plus),
        'b_enclosure': str(b),
        'gamma_delta_degrees': [str(gamma*180/pi), str(delta*180/pi)],
        'area_ball': str(area),
        'strict_global_optimum_bracket': [area_lower, area_upper],
        'scope': 'Global minimum over both arbitrary translations at theta_5=0 only.'
    }, indent=2))


if __name__ == '__main__':
    main()
