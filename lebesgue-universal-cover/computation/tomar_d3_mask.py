"""Exact raw-node mask for the triangular dihedral normalization."""
from fractions import Fraction
import numpy as np

CUT = Fraction(20943951024, 10**11)


def leaf_mask(codes, part, angle_high):
    """Values 1, 2, 3 denote inside, crossing, strictly outside. Internal nodes are zero."""
    assert codes.dtype == np.uint8 and codes.ndim == 1
    # Scale by the numerator of angle_high and CUT's denominator. The upper
    # comparison endpoint doubles on each split in the angle coordinate.
    stride = angle_high.numerator * CUT.denominator
    threshold = 2 * CUT.numerator * angle_high.denominator
    root_bit = (part >> 8) & 1  # bit for coordinate 2, C-order over eleven axes
    stack = [(root_bit * stride, (root_bit + 1) * stride, threshold)]
    mask = np.zeros(len(codes), dtype=np.uint8)
    for pos, code in enumerate(codes.tolist()):
        lo, hi, bound = stack.pop()
        if code >= 253:
            mask[pos] = 1 if hi <= bound else 3 if lo > bound else 2
        elif code == 2:
            mid = lo + hi
            stack.append((mid, 2 * hi, 2 * bound))
            stack.append((2 * lo, mid, 2 * bound))
        else:
            assert 0 <= code < 11
            stack.extend(((lo, hi, bound), (lo, hi, bound)))
    assert not stack
    return mask


