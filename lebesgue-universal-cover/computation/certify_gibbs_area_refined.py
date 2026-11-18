"""Outward rational Green integral of the Gibbs R34 boundary at one rational angle."""
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys

import rational_trig as trig

ROUND = 10**48
trig.T = ROUND
pi_machin = trig.add(trig.scale(trig.arctan_inv(5, 64), F(16)),
                     trig.scale(trig.arctan_inv(239, 20), F(-4)))


class I:
    def __init__(self, lo, hi=None):
        if isinstance(lo, I):
            self.lo, self.hi = lo.lo, lo.hi
        else:
            self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def round(lo, hi):
        return I(F((lo * ROUND).__floor__(), ROUND), F((hi * ROUND).__ceil__(), ROUND))

    def __add__(self, other):
        b = I(other)
        return I.round(self.lo + b.lo, self.hi + b.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -I(other)

    def __rsub__(self, other):
        return I(other) + -self

    def __mul__(self, other):
        b = I(other)
        values = [u * v for u in (self.lo, self.hi) for v in (b.lo, b.hi)]
        return I.round(min(values), max(values))

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = I(other)
        assert b.lo * b.hi > 0, b.out()
        return self * I.round(1 / b.hi, 1 / b.lo)

    def __rtruediv__(self, other):
        return I(other) / self

    def sqrt(self):
        assert self.lo >= 0, self.out()
        lo = F(isqrt((self.lo * ROUND**2).__floor__()), ROUND)
        hi = F(isqrt((self.hi * ROUND**2).__floor__()) + 1, ROUND)
        assert lo * lo <= self.lo and hi * hi >= self.hi
        return I(lo, hi)

    def out(self):
        return [str(self.lo), str(self.hi)]


if len(sys.argv) not in (2, 3, 4):
    raise SystemExit("usage: certify_gibbs_area_refined.py SIGMA_DEGREES [REMAINDER_ORDER=8] [CELLS=128]")
sigma = F(sys.argv[1])
remainder_order = int(sys.argv[2]) if len(sys.argv) > 2 else 8
cells = int(sys.argv[3]) if len(sys.argv) > 3 else 128
assert 0 < sigma < 30 and remainder_order >= 2 and remainder_order % 2 == 0 and cells >= 1
order = remainder_order + 1


class Jet:
    """Ordinary Taylor coefficients about an arbitrary centre in the input interval."""
    def __init__(self, value):
        if isinstance(value, Jet):
            self.c = value.c[:]
        else:
            self.c = [I(value)] + [I(0) for _ in range(order)]

    @staticmethod
    def from_coefficients(values):
        result = Jet(0)
        result.c = values
        return result

    def __add__(self, other):
        b = Jet(other)
        return Jet.from_coefficients([a + z for a, z in zip(self.c, b.c)])

    __radd__ = __add__

    def __neg__(self):
        return Jet.from_coefficients([-a for a in self.c])

    def __sub__(self, other):
        return self + -Jet(other)

    def __rsub__(self, other):
        return Jet(other) + -self

    def __mul__(self, other):
        b = Jet(other)
        return Jet.from_coefficients([
            sum(self.c[k] * b.c[j-k] for k in range(j+1)) for j in range(order+1)
        ])

    __rmul__ = __mul__

    def reciprocal(self):
        a = self.c
        b = [1 / a[0]]
        for j in range(1, order+1):
            b.append(-sum(a[k] * b[j-k] for k in range(1, j+1)) / a[0])
        return Jet.from_coefficients(b)

    def __truediv__(self, other):
        return self * Jet(other).reciprocal()

    def __rtruediv__(self, other):
        return Jet(other) * self.reciprocal()

    def sqrt(self):
        a = self.c
        b = [a[0].sqrt()]
        for j in range(1, order+1):
            b.append((a[j] - sum(b[k] * b[j-k] for k in range(1, j))) / (2*b[0]))
        return Jet.from_coefficients(b)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def mul(v, a):
    return tuple(v * x for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def perp(a):
    return (-a[1], a[0])


half = F(1, 2)
r3 = I(3).sqrt()


def meet(a, b, ra=half, rb=half):
    d = cross(a, b)
    return ((ra*b[1] - rb*a[1])/d, (rb*a[0] - ra*b[0])/d)


def project(x, n):
    return sub(x, mul(dot(x, n)-half, n))


def cl(c, n, sign):
    d = dot(n, c) - half
    return add(project(c, n), mul(sign * (1-d*d).sqrt(), perp(n)))


def cc(a, b, sign):
    q = sub(b, a)
    return add(mul(half, add(a, b)), mul(sign * (1/dot(q,q)-F(1,4)).sqrt(), perp(q)))


def atan(z):
    bound = max(abs(z.lo), abs(z.hi))
    assert bound < half
    result = I(0)
    power = z
    terms = 80
    for k in range(terms):
        result = result + (-1)**k * power/(2*k+1)
        power = power * z * z
    omitted = bound**(2*terms+1) / (2*terms+1)
    return result + I(-omitted, omitted)


def arc(a, b, c):
    u, v = sub(a,c), sub(b,c)
    denominator = dot(u,v)
    assert denominator.lo > 0 and cross(u,v).lo > 0
    return cross(c, sub(b,a)) + atan(cross(u,v)/denominator)


n = [(I(1), I(0)), (I(half),r3/2), (I(-half),r3/2),
     (I(-1),I(0)), (I(-half),-r3/2), (I(half),-r3/2)]
angle_m0 = trig.scale(pi_machin, F(1,6) + sigma/180)
m0 = I(*trig.trig(angle_m0, "cos")), I(*trig.trig(angle_m0, "sin"))
m = [(a*m0[0]-b*m0[1], b*m0[0]+a*m0[1]) for a,b in n]
b1, f1, d1 = meet(n[0], n[1]), meet(n[2], n[3]), meet(n[4], n[5])
b3 = meet(n[1], m[0])
c3 = meet(n[0], m[5])
e3 = meet(n[4], m[3])
f3 = meet(n[3], m[2])
g = cl(f3, n[5], 1)
q, p, ip, h, k = (project(e3,n[1]), project(g,n[2]), project(b3,m[3]),
                   project(c3,n[3]), project(f3,m[5]))
r = meet(m[2], n[1], half, 0)
s = cl(r,n[0],-1)
t_eh = cc(s,b3,1)
u = cl(c3,(-r3/2,I(-half)),-1)

c = (-e3[0],e3[1])
mid = (I(F(-1,4)),-r3/4)
l0 = (r3/2-half,(r3/2-half-1)/r3)
z = cc(l0,mid,-1)
y = (I(0),c[1]+(1-c[0]*c[0]).sqrt())
v = sub(z,c)
v2 = dot(v,v)
t = add(c,sub(mul(1/v2,v),mul((v2-1).sqrt()/v2,perp(v))))
direction = sub(t,z)
n_l, n_n = (half,-r3/2), (-half,-r3/2)
e_l, e_n = (r3/2,half), (-r3/2,half)


def pair(point):
    cast = Jet if isinstance(point[0], Jet) else I
    make_vec = lambda arr: tuple(cast(value) for value in arr)

    def inverse(normal, along):
        normal, along = make_vec(normal), make_vec(along)
        d = dot(normal,point)-half
        return sub(sub(point,mul(d,normal)),mul((1-d*d).sqrt(),along))

    l = inverse(n_l,e_l)
    nn = inverse(n_n,e_n)
    v0 = l[0]+half
    lx = (v0 + cast(r3)*(1-v0*v0).sqrt())/2-half
    return cc((lx,(lx-1)/cast(r3)),sub(mul(2,make_vec(mid)),nn),-1)


def density_coefficients(position):
    w = Jet(position)
    w.c[1] = I(1)
    x, y = pair(add(tuple(map(Jet,z)),mul(w,tuple(map(Jet,direction)))))
    return [
        sum((j+1)*(x.c[k-j]*y.c[j+1]-y.c[k-j]*x.c[j+1]) for j in range(k+1))
        for k in range(remainder_order+1)
    ]


width = F(1,cells)
denominators = {k: (k+1)*2**k for k in range(0,remainder_order+1,2)}
factors = {k: width**(k+1)/d for k,d in denominators.items()}
integral = I(0)
maximum_remainder_coefficient = F(0)
for j in range(cells):
    interval = I(j*width,(j+1)*width)
    centre = I((2*j+1)*width/2)
    f_centre = density_coefficients(centre)
    f_remainder = density_coefficients(interval)[remainder_order]
    integral += sum(factors[k]*f_centre[k] for k in factors if k < remainder_order)
    integral += factors[remainder_order]*f_remainder
    maximum_remainder_coefficient = max(maximum_remainder_coefficient,
                                        abs(f_remainder.lo),abs(f_remainder.hi))

tp = pair(t)
fixed = (cross(b1,q)+arc(q,y,e3)+arc(y,t,c)+cross(t,z)
         +arc(tp,p,g)+cross(p,f1)+cross(f1,h)+arc(h,u,c3)
         +cross(u,t_eh)+arc(t_eh,ip,b3)+cross(ip,e3)+cross(e3,d1)
         +cross(d1,g)+arc(g,k,f3)+cross(k,c3)+cross(c3,b1))
area = (fixed+integral)/2
assert area.lo > 0
print(json.dumps({
    "python": sys.version,
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "trig_sha256": hashlib.sha256(Path(trig.__file__).read_bytes()).hexdigest(),
    "sigma_degrees": str(sigma),
    "rounding_denominator": str(ROUND),
    "cells": cells,
    "Taylor_remainder_order": remainder_order,
    "maximum_remainder_coefficient": str(maximum_remainder_coefficient),
    "integral_remainder_abs_upper": str(cells*factors[remainder_order]*maximum_remainder_coefficient),
    "quadrature_interval": integral.out(),
    "fixed_Green_interval": fixed.out(),
    "area_interval": area.out(),
    "area_width": str(area.hi-area.lo),
    "passed": True,
},indent=2))
