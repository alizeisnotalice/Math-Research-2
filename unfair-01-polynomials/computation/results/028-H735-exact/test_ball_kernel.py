"""Corner-enumeration tests for directed centre-radius multiplication."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib, json
from verify_h735_ball_certificate import balls_from_boxes, coefficient, multiply


def exact_product(polynomials):
    answer = [F(1)]
    for polynomial in polynomials:
        next_answer = [F(0)] * (len(answer) + len(polynomial) - 1)
        for i, a in enumerate(answer):
            for j, b in enumerate(polynomial):
                next_answer[i + j] += a * b
        answer = next_answer
    return answer


def test_case(intervals, bits, cutoff):
    D, factors = balls_from_boxes({i: p for i, p in enumerate(intervals)}, bits)
    result = [(D, 0)]
    for i in range(len(intervals)):
        result = multiply(result, factors[i], D, cutoff)
    coordinate_intervals = [item for polynomial in intervals for item in polynomial]
    for choices in product((0, 1), repeat=len(coordinate_intervals)):
        it = iter(choices)
        polynomials = [[pair[next(it)] for pair in polynomial] for polynomial in intervals]
        exact = exact_product(polynomials)
        for degree, value in enumerate(exact[:cutoff + 1]):
            centre, radius = result[degree]
            assert F(centre - radius, D) <= value <= F(centre + radius, D)


def boundary_rounding_test():
    """A sub-grid negative singleton must not acquire a false strict sign."""
    D, factors = balls_from_boxes({0: [(F(0), F(0))]}, 96)
    centre, radius = coefficient(factors, [0], 0, D)
    assert (centre, radius) == (0, 0)

    # -2^-100 lies inside one 96-bit grid unit of zero.  Its outward ball
    # has upper endpoint zero, hence the strict-negative certificate rejects.
    D, factors = balls_from_boxes({0: [(F(-1, 1 << 100), F(-1, 1 << 100))]}, 96)
    centre, radius = coefficient(factors, [0], 0, D)
    assert centre + radius == 0

    # -2^-90 is exactly 64 grid units below zero and is strictly certified.
    D, factors = balls_from_boxes({0: [(F(-1, 1 << 90), F(-1, 1 << 90))]}, 96)
    centre, radius = coefficient(factors, [0], 0, D)
    assert centre + radius < 0


def main():
    # Both signs, zero-crossing, and a negative nonintegral centre exercise
    # min/max products and Python floor division in the round-down branch.
    cases = [
        ([[(F(-5, 7), F(2, 3)), (F(-3, 2), F(-1, 4)), (F(1), F(1))],
          [(F(-2), F(-1)), (F(0), F(5, 6)), (F(1), F(1))]], 7, 4),
        ([[(F(-1, 9), F(1, 7)), (F(-11, 13), F(-3, 11))],
          [(F(-2, 5), F(4, 9)), (F(-7, 8), F(-1, 6))],
          [(F(1, 3), F(5, 4)), (F(-1, 2), F(2, 3))]], 5, 2),
    ]
    for args in cases:
        test_case(*args)
    boundary_rounding_test()
    result = {'status': 'corner_enumeration_passed', 'cases': len(cases),
              'boundary_cases': ['zero', '-2^-100_not_strict', '-2^-90_strict'],
              'notes': 'All endpoint combinations, including sign-crossing and cutoff terms, lie in directed balls.'}
    Path(__file__).with_name('ball-kernel-tests.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
