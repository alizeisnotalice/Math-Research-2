from math import isqrt
def req(x,s):
 if not x: raise ValueError(s)

def _trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def _mul_mod(left, right, prime):
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] = (out[i + j] + a * b) % prime
    return _trim(out)


def _rem_mod(left, right, prime):
    left, right = _trim([x % prime for x in left]), _trim([x % prime for x in right])
    req(right != [0], 'finite-field division by zero')
    inverse = pow(right[-1], -1, prime)
    while left != [0] and len(left) >= len(right):
        shift = len(left) - len(right)
        coefficient = left[-1] * inverse % prime
        for j, value in enumerate(right):
            left[shift + j] = (left[shift + j] - coefficient * value) % prime
        left = _trim(left)
    return left


def _pow_mod(base, exponent, modulus, prime):
    result = [1]
    base = _rem_mod(base, modulus, prime)
    while exponent:
        if exponent & 1:
            result = _rem_mod(_mul_mod(result, base, prime), modulus, prime)
        base = _rem_mod(_mul_mod(base, base, prime), modulus, prime)
        exponent //= 2
    return result


def _gcd_mod(left, right, prime):
    while right != [0]:
        left, right = right, _rem_mod(left, right, prime)
    inverse = pow(left[-1], -1, prime)
    return [(x * inverse) % prime for x in left]


def _prime(value):
    return type(value) is int and value >= 2 and all(value % divisor for divisor in range(2, isqrt(value) + 1))


def _irreducible_mod(poly, prime):
    """Rabin's exact irreducibility test over F_prime."""
    degree = len(poly) - 1
    req(degree >= 1 and poly[-1] == 1 and all(type(x) is int and 0 <= x < prime for x in poly),
        'monic finite-field factor encoding')
    divisors = [q for q in range(2, degree + 1)
                if degree % q == 0 and _prime(q)]
    x = _rem_mod([0, 1], poly, prime)
    frobenius = x
    for step in range(1, degree + 1):
        frobenius = _pow_mod(frobenius, prime, poly, prime)
        if step in {degree // q for q in divisors}:
            difference = _trim([(frobenius[i] if i < len(frobenius) else 0) -
                                (x[i] if i < len(x) else 0)
                                for i in range(max(len(frobenius), len(x)))])
            req(_gcd_mod(poly, difference, prime) == [1], 'claimed modular factor reducible')
    req(frobenius == x, 'Rabin final Frobenius test')
