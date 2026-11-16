"""Search-independent 96-bit directed balls for top-coefficient witnesses."""
from fractions import Fraction as F


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


def balls_from_boxes(boxes, bits=96):
    """Convert ascending Fraction coefficient boxes into integer D-balls."""
    scale = 1 << bits
    output = {}
    for index, polynomial in boxes.items():
        row = []
        for lower, upper in polynomial:
            assert lower <= upper
            midpoint, halfwidth = (lower + upper) / 2, (upper - lower) / 2
            centre = (midpoint.numerator * scale) // midpoint.denominator
            radius = ceil_fraction((halfwidth + abs(midpoint - F(centre, scale))) * scale)
            assert F(centre - radius, scale) <= lower <= upper <= F(centre + radius, scale)
            row.append((centre, radius))
        output[index] = row
    return scale, output


def multiply(left, right, scale, cutoff):
    """Directed product truncated only beyond the requested reversed degree."""
    length = min(cutoff + 1, len(left) + len(right) - 1)
    centres, errors = [0] * length, [0] * length
    for i, (centre, radius) in enumerate(left):
        for j, (other, other_radius) in enumerate(right):
            if i + j >= length:
                break
            centres[i + j] += centre * other
            errors[i + j] += abs(centre) * other_radius + abs(other) * radius + radius * other_radius
    result = []
    for unrounded, error in zip(centres, errors):
        centre = unrounded // scale
        rounding_remainder = unrounded - centre * scale
        radius = (error + rounding_remainder + scale - 1) // scale
        result.append((centre, radius))
    return result


def top_coefficient(factors, indices, top_offset, scale):
    """Certify coefficient degree(total_degree-top_offset) by reversing boxes."""
    total_degree = sum(len(factors[index]) - 1 for index in indices)
    assert type(top_offset) is int and 0 <= top_offset <= total_degree
    result = [(scale, 0)]
    for index in indices:
        result = multiply(result, list(reversed(factors[index])), scale, top_offset)
    return total_degree - top_offset, result[top_offset]
