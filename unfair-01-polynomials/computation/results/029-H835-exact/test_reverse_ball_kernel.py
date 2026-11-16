"""Endpoint and reversal tests for H8 top-offset directed-ball arithmetic."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

from reverse_ball_kernel import balls_from_boxes, multiply, top_coefficient


def exact_product(polynomials):
    result = [F(1)]
    for polynomial in polynomials:
        next_result = [F(0)] * (len(result) + len(polynomial) - 1)
        for i, left in enumerate(result):
            for j, right in enumerate(polynomial):
                next_result[i + j] += left * right
        result = next_result
    return result


def main():
    intervals = [
        [(F(-3, 7), F(2, 5)), (F(-5, 6), F(-1, 8)), (F(1), F(1))],
        [(F(-2, 3), F(1, 9)), (F(-7, 5), F(-3, 10)), (F(1), F(1))],
        [(F(-1, 4), F(3, 8)), (F(1, 3), F(4, 5)), (F(1), F(1))],
    ]
    scale, balls = balls_from_boxes({index: polynomial for index, polynomial in enumerate(intervals)}, 11)
    reversed_product = [(scale, 0)]
    for index in range(len(intervals)):
        reversed_product = multiply(reversed_product, list(reversed(balls[index])), scale, 3)

    coordinates = [interval for polynomial in intervals for interval in polynomial]
    for choices in product((0, 1), repeat=len(coordinates)):
        iterator = iter(choices)
        polynomials = [[interval[next(iterator)] for interval in polynomial] for polynomial in intervals]
        exact = exact_product(polynomials)
        for offset in range(4):
            degree, (centre, radius) = top_coefficient(balls, range(3), offset, scale)
            assert degree == len(exact) - 1 - offset
            assert F(centre - radius, scale) <= exact[degree] <= F(centre + radius, scale)

    # A monic quadratic has reversed leading coefficient one.  This separately
    # locks the convention used by compact top-offset witnesses.
    assert list(reversed(balls[0]))[0] == (scale, 0)
    result = {'status': 'reverse_directed_ball_kernel_passed',
              'endpoint_combinations': 2 ** len(coordinates),
              'top_offsets_checked': 4,
              'notes': 'Reversed boxes certify original degree total_degree-offset; truncation is only above that reversed degree.'}
    Path(__file__).with_name('reverse-ball-kernel-tests.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
