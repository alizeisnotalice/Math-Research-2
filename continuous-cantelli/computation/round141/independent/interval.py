"""Directed Decimal bounds for beta-weight heat densities on spatial boxes."""
import argparse
import json
from decimal import Context, Decimal as D, ROUND_CEILING, ROUND_FLOOR
from pathlib import Path

LO = Context(prec=48, rounding=ROUND_FLOOR)
HI = Context(prec=48, rounding=ROUND_CEILING)
NEAR = Context(prec=48)


class Box:
    def __init__(self, a, b=None):
        self.a = D(a)
        self.b = D(a if b is None else b)

    def __add__(self, other):
        return Box(LO.add(self.a, other.a), HI.add(self.b, other.b))

    def __sub__(self, other):
        return Box(LO.subtract(self.a, other.b), HI.subtract(self.b, other.a))

    def __mul__(self, other):
        pairs = [(x, y) for x in (self.a, self.b) for y in (other.a, other.b)]
        return Box(min(LO.multiply(x, y) for x, y in pairs),
                   max(HI.multiply(x, y) for x, y in pairs))

    def __truediv__(self, other):
        assert other.a > 0
        return self * Box(LO.divide(1, other.b), HI.divide(1, other.a))

    def exp(self):
        return Box(NEAR.next_minus(NEAR.exp(self.a)), NEAR.next_plus(NEAR.exp(self.b)))

    def sqrt(self):
        return Box(NEAR.next_minus(NEAR.sqrt(self.a)), NEAR.next_plus(NEAR.sqrt(self.b)))

    def pow(self, n):
        assert 0 <= self.a
        m = Box(1)
        p = self
        while n:
            if n % 2:
                m = m * p
            n //= 2
            if n:
                p = p * p
        return m

    def json(self):
        return [str(self.a), str(self.b)]


ZERO, ONE, TWO = Box(0), Box(1), Box(2)
SHIFT0, SHIFT1 = D('0.375'), D('0.875')
E = Box('0.6') * ((ZERO - Box(16)) / Box(3)).exp()


def gaussian(v, y):
    if v == 0:
        assert y > 0
        return ZERO
    return ((ZERO - Box(y) * Box(y)) / (TWO * Box(v))).exp() / Box(v).sqrt()


def weight(t, j, k):
    if t == 0:
        return ZERO
    a = (ZERO - Box(t)).exp()
    return (ZERO - Box(j) * Box(t)).exp() * (ONE - a).pow(k)


def gaussian_box(lo, hi, left, right):
    # For y>0 the Gaussian decreases with y. In variance it has just one maximum.
    assert 0 <= lo < hi and 0 < left <= right
    critical = left*left
    lower = min(gaussian(lo, right).a, gaussian(hi, right).a)
    upper = max(gaussian(lo, left).b, gaussian(hi, left).b)
    if lo < critical < hi:
        upper = gaussian(critical, left).b
    return Box(lower, upper)


def certify(j, k, left, right, depth, model='three'):
    assert 0 < j and 0 < k and D(0) < left < right
    if model == 'three':
        assert D(2) < left
    center = D(2) if model == 'three' else D(0)
    mass = E if model == 'three' else Box('0.55')
    critical = Box(j+k) / Box(j)
    w_peak = (Box(j) / Box(j+k)).pow(j) * (Box(k) / Box(j+k)).pow(k)
    width = D(1) / (2 ** (depth+1))
    results = []
    for segment in range(3 if model == 'three' else 2):
        bound = ZERO
        for i in range(2**depth):
            lo = D(segment)/2 + i*width
            hi = lo + width
            wl = weight(lo, j, k)
            wh = weight(hi, j, k)
            if Box(hi).exp().b <= critical.a or Box(lo).exp().a >= critical.b:
                w = Box(min(wl.a, wh.a), max(wl.b, wh.b))
            else:
                w = Box(min(wl.a, wh.a), w_peak.b)
            if segment != 2:
                free = gaussian_box(lo, hi, left, right)
            if segment == 0:
                q = free
            elif segment == 1:
                source = gaussian_box(lo-SHIFT0, hi-SHIFT0, left-center, right-center)
                q = free - mass * source
            else:
                source = gaussian_box(lo-SHIFT1, hi-SHIFT1, left-center, right-center)
                q = mass * source
            bound = bound + w * q * Box(width)
        results.append(bound)
    total = sum(results, ZERO)
    payload = dict(model=model, j=j, k=k, spatial_interval=[str(left), str(right)], cells_per_segment=2**depth,
        # The enclosed quantity is sqrt(2*pi) times the density of D.
        certified_stripped_density=total.json(), normalized_by_peak=(total/w_peak).json(),
        stripped_segment_bounds=[x.json() for x in results], peak=w_peak.json())
    return total, payload


def main():
    p = argparse.ArgumentParser()
    p.add_argument('j', type=int)
    p.add_argument('k', type=int)
    p.add_argument('left', type=D)
    p.add_argument('right', type=D)
    p.add_argument('depth', type=int)
    p.add_argument('out', type=Path)
    p.add_argument('--model', choices=('two','three'), default='three')
    args = p.parse_args()
    box, record = certify(args.j, args.k, args.left, args.right, args.depth, args.model)
    args.out.write_text(json.dumps(record, indent=2) + '\n')
    print(record['normalized_by_peak'])
    assert box.b < 0


if __name__ == '__main__':
    main()
