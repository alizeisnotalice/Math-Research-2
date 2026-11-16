from fractions import Fraction as F
from math import comb, isqrt


def pf(value):
    numerator, denominator = value.split('/', 1)
    return F(int(numerator), int(denominator))


def cj(value):
    return pf(value['re']), pf(value['im'])


def sub(left, right):
    return left[0] - right[0], left[1] - right[1]


def mul(left, right):
    return left[0] * right[0] - left[1] * right[1], left[0] * right[1] + left[1] * right[0]


def pw(value, exponent):
    result = (F(1), F(0))
    while exponent:
        if exponent & 1:
            result = mul(result, value)
        value = mul(value, value)
        exponent //= 2
    return result


def ab(value, scale):
    squared = value[0] * value[0] + value[1] * value[1]
    lower_square = (squared.numerator * scale * scale) // squared.denominator
    lower = isqrt(lower_square)
    return F(lower, scale), F(lower + 1, scale)


def tay(coefficients, point):
    output = []
    for derivative in range(len(coefficients)):
        value = (F(0), F(0))
        for degree in range(derivative, len(coefficients)):
            power = pw(point, degree - derivative)
            value = (value[0] + coefficients[degree] * comb(degree, derivative) * power[0],
                     value[1] + coefficients[degree] * comb(degree, derivative) * power[1])
        output.append(value)
    return output


def mother(coefficients, point, radius, scale):
    taylor = tay(coefficients, point)
    _, upper_constant = ab(taylor[0], scale)
    for degree in range(2, len(taylor)):
        _, upper = ab(taylor[degree], scale)
        upper_constant += upper * radius ** degree
    lower_linear, _ = ab(taylor[1], scale)
    return upper_constant, lower_linear * radius


def lift(parent, parent_radius, point, radius, scale, exponent):
    _, upper_constant = ab(sub(pw(point, exponent), parent), scale)
    upper_constant += parent_radius
    for degree in range(2, exponent + 1):
        _, upper = ab((comb(exponent, degree) * pw(point, exponent - degree)[0],
                       comb(exponent, degree) * pw(point, exponent - degree)[1]), scale)
        upper_constant += upper * radius ** degree
    lower_linear, _ = ab((exponent * pw(point, exponent - 1)[0],
                           exponent * pw(point, exponent - 1)[1]), scale)
    return upper_constant, lower_linear * radius


def dis(left, left_radius, right, right_radius):
    return (left[0] - right[0]) ** 2 + (left[1] - right[1]) ** 2 > (left_radius + right_radius) ** 2
