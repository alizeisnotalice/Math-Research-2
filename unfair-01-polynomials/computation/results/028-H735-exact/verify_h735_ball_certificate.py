"""Independent directed centre-radius verifier for compact H7(3,5) data.

No search implementation or search coefficient boxes are imported.  Root boxes
are rebuilt by the copied independent front end, then converted to 96-bit
integer balls.  A compact certificate supplies only a side and coefficient
witness for each necessary root mask.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, sys, time
from math import isqrt
from verify_frontend import verify, H, sha, ALL, req


def ceil_fraction(x):
    return -((-x.numerator) // x.denominator)


# This deliberately duplicates only the elementary finite-field checks needed
# here.  It does not import the searcher, SymPy, or a prior "passed" receipt.
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


def verify_modular_galois(galois_path):
    """Independently check the fixed H and its S7 modular certificate.

    The three square-free reductions certify a 7-cycle, a 6-cycle, and a
    (5)(2) cycle.  The standard group argument used by the project turns
    these into S7; this checker establishes the integer premises only.
    """
    raw = json.loads(galois_path.read_text())
    req(isinstance(raw, dict) and set(raw) == {'certificates'}, 'galois outer schema')
    certificates = raw['certificates']
    req(isinstance(certificates, list) and len(certificates) == 1, 'one fixed H galois certificate')
    item = certificates[0]
    req(type(item.get('degree')) is int and item['degree'] == 7, 'H degree')
    req(item.get('coefficients_ascending') == H and
        all(type(value) is int for value in item['coefficients_ascending']), 'galois certificate fixed H')
    req(H != list(reversed(H)), 'fixed H must be nonreciprocal')
    req(item.get('name') == 'H7_1_x_x2_x6_x7', 'fixed H name')
    reductions = item.get('reductions')
    req(isinstance(reductions, list) and len(reductions) == 3, 'three modular reductions')
    observed = {}
    for row in reductions:
        req(isinstance(row, dict) and set(row) == {'prime', 'factors'}, 'modular reduction schema')
        prime = row['prime']
        req(_prime(prime) and prime not in observed, 'unique prime reduction')
        factors = row['factors']
        req(isinstance(factors, list) and factors, 'nonempty modular factorization')
        product = [1]
        degrees = []
        for factor in factors:
            req(isinstance(factor, dict) and set(factor) == {'coefficients_ascending', 'multiplicity'},
                'modular factor schema')
            req(type(factor['multiplicity']) is int and factor['multiplicity'] == 1,
                'squarefree claimed factor multiplicity')
            polynomial = factor['coefficients_ascending']
            _irreducible_mod(polynomial, prime)
            product = _mul_mod(product, polynomial, prime)
            degrees.append(len(polynomial) - 1)
        req(product == _trim([value % prime for value in H]), 'modular factor product')
        derivative = _trim([(i * H[i]) % prime for i in range(1, len(H))])
        req(_gcd_mod(product, derivative, prime) == [1], 'modular reduction not squarefree')
        observed[prime] = sorted(degrees)
    req(observed == {7: [7], 5: [1, 6], 11: [2, 5]}, 'fixed S7 cycle types')
    return {'status': 'fixed_H_S7_modular_premises_verified_by_integer_arithmetic',
            'galois_sha256': sha(galois_path), 'degree': 7,
            'nonreciprocal': True, 'reduction_degrees': {str(key): value for key, value in observed.items()}}


def balls_from_boxes(boxes, bits):
    """Return balls (integer centre, integer radius) at scale D=2**bits.

    For [lo,hi], m=(lo+hi)/2 and rho=(hi-lo)/2, c=floor(D*m),
    r=ceil(D*(rho+abs(m-c/D))).  Thus [lo,hi] is contained in
    [(c-r)/D,(c+r)/D], even with a negative centre.
    """
    D = 1 << bits
    out = {}
    for index, polynomial in boxes.items():
        row = []
        for lo, hi in polynomial:
            assert lo <= hi
            mid, radius = (lo + hi) / 2, (hi - lo) / 2
            centre = (mid.numerator * D) // mid.denominator
            error = radius + abs(mid - F(centre, D))
            bound = ceil_fraction(error * D)
            assert F(centre - bound, D) <= lo <= hi <= F(centre + bound, D)
            row.append((centre, bound))
        out[index] = row
    return D, out


def multiply(left, right, D, cutoff):
    """Ball product, retaining exactly coefficients at most cutoff.

    If x=C±R and y=M±S at scale D, then xy has centre CM/D and
    error (|C|S+|M|R+RS)/D.  We replace its centre by floor(CM/D), adding
    the corresponding downward-rounding remainder to its radius.
    """
    length = min(cutoff + 1, len(left) + len(right) - 1)
    centre = [0] * length
    error = [0] * length
    for i, (c, r) in enumerate(left):
        for j, (m, s) in enumerate(right):
            if i + j >= length:
                break
            centre[i + j] += c * m
            error[i + j] += abs(c) * s + abs(m) * r + r * s
    result = []
    for value, uncertainty in zip(centre, error):
        rounded = value // D
        remainder = value - rounded * D
        radius = (uncertainty + remainder + D - 1) // D
        result.append((rounded, radius))
    return result


def coefficient(factors, indices, degree, D):
    result = [(D, 0)]
    for index in indices:
        result = multiply(result, factors[index], D, degree)
    return result[degree]


def fraction_coefficient(factors, indices, degree):
    """Raw Fraction interval convolution for selected final witness samples."""
    lower, upper = [F(1)], [F(1)]
    for index in indices:
        factor = factors[index]
        length = min(degree + 1, len(lower) + len(factor) - 1)
        next_lower, next_upper = [F(0)] * length, [F(0)] * length
        for i, (left, right) in enumerate(zip(lower, upper)):
            for j, (lo, hi) in enumerate(factor):
                if i + j >= length:
                    break
                products = (left * lo, left * hi, right * lo, right * hi)
                next_lower[i + j] += min(products)
                next_upper[i + j] += max(products)
        lower, upper = next_lower, next_upper
    return lower[degree], upper[degree]


def standard_mask(balls):
    return sum(1 << x['index'] for x in balls['groups'] if x['family'] == 'a')


def compact_records(data, layouts, standard):
    req(data['format'] == 'H7-a3b5-compact-ball-witness-v1', 'certificate format')
    req(data['status'] == 'complete', 'certificate complete status')
    witnesses = data['witnesses']
    req(isinstance(witnesses, list) and len(witnesses) == len(layouts) == 119062, 'witness count')
    out = []
    for (energy, mask), item in zip(layouts, witnesses):
        req(isinstance(item, list) and len(item) == 2, 'witness row shape')
        side, degree = item
        req(type(side) is int and type(degree) is int, 'witness integer encoding')
        if mask == standard:
            req(energy == 15 and side == degree == -1, 'unique exact standard witness')
        else:
            req(side in (0, 1) and degree >= 0, 'nonstandard witness encoding')
        out.append((energy, mask, side, degree))
    return out


def validate_certificate_metadata(certificate, balls_path, canonical_sha256):
    """Reject a compact stream attached to the wrong roots/precision."""
    req(certificate['balls_sha256'] == sha(balls_path), 'balls sha binding')
    req(certificate['canonical_sha256'] == canonical_sha256, 'canonical sha binding')
    # Compact data records the 128-bit search precision.  This verifier uses
    # its own coarser 96-bit balls below, deliberately independent of that
    # search arithmetic.
    req(type(certificate['bits']) is int and certificate['bits'] == 128,
        'expected integer 128-bit source certificate')


def check(balls_path, masks_path, certificate_path, galois_path, record_path):
    balls = json.loads(balls_path.read_text())
    front, layouts, boxes = verify(balls)
    galois = verify_modular_galois(galois_path)
    certificate = json.loads(certificate_path.read_text())
    validate_certificate_metadata(certificate, balls_path, front['canonical_sha256'])
    expected = set(layouts)
    external = [tuple(x) for x in json.loads(masks_path.read_text())]
    req(set(external) == expected and len(external) == len(expected), 'external canonical mask set')
    standard = standard_mask(balls)
    standard_indices = [i for i in range(29) if standard >> i & 1]
    complement_indices = [i for i in range(29) if not (standard >> i & 1)]
    # These are root degrees, not coefficient indices.  They separately bind
    # the unique binary a/b split used for the standard row.
    req(sum(len(boxes[i]) - 1 for i in standard_indices) == 21,
        'standard a factor root degree 21')
    req(sum(len(boxes[i]) - 1 for i in complement_indices) == 35,
        'standard complement root degree 35')
    entries = compact_records(certificate, layouts, standard)
    got = {(energy, mask) for energy, mask, _, _ in entries}
    req(got == expected and len(got) == len(entries), 'certificate exact complete coverage')
    D, factors = balls_from_boxes(boxes, 96)
    started = time.monotonic()
    minimum = None
    checked_negative = 0
    raw_sample_keys, raw_samples = set(), []
    for n, (energy, mask, side, degree) in enumerate(entries, 1):
        if mask == standard:
            continue
        req(side in (0, 1), 'nonstandard side')
        selected = mask if side == 0 else ALL ^ mask
        indices = [i for i in range(29) if selected >> i & 1]
        max_degree = sum(len(factors[i]) - 1 for i in indices)
        req(0 <= degree <= max_degree, 'witness coefficient degree')
        centre, radius = coefficient(factors, indices, degree, D)
        req(centre + radius < 0, 'independent ball coefficient not strictly negative')
        sample_key = (energy, side)
        if sample_key not in raw_sample_keys and len(raw_samples) < 12:
            raw_sample_keys.add(sample_key)
            raw_lo, raw_hi = fraction_coefficient(boxes, indices, degree)
            req(raw_hi < 0, 'raw Fraction sample not strictly negative')
            raw_samples.append({'E': energy, 'side': side, 'mask': mask, 'index': degree,
                                'raw_lo': str(raw_lo), 'raw_hi': str(raw_hi)})
        # floor((-upper)*D)/D, a conservative lower bound for -upper.
        margin = -(centre + radius)
        req(margin > 0, 'positive integer margin')
        minimum = margin if minimum is None else min(minimum, margin)
        checked_negative += 1
        if n % 1000 == 0:
            record_path.with_suffix('.progress.json').write_text(json.dumps(
                {'status': 'verifying', 'checked_rows': n,
                 'negative_rows': checked_negative,
                 'elapsed_seconds': time.monotonic() - started}) + '\n')
    req(checked_negative == 119061 and minimum is not None, 'negative/standard partition')
    result = {
        'status': 'all_119062_compact_witnesses_verified_by_independent_balls',
        'negative_layouts': checked_negative,
        'standard_layouts': 1,
        'standard_mask': standard,
        'standard_factor_root_degrees': {'a': 21, 'b': 35},
        'minimum_certified_negative_margin': str(F(minimum, D)),
        'raw_fraction_samples': raw_samples,
        'bits': 96,
        'arithmetic': 'independent integer centre-radius balls from Fraction root boxes; cutoff only above requested degree',
        'elapsed_seconds': time.monotonic() - started,
        'balls_sha256': sha(balls_path),
        'canonical_masks_sha256': sha(masks_path),
        'certificate_sha256': sha(certificate_path),
        'galois': galois,
        'frontend_sha256': sha(Path(__file__).with_name('verify_frontend.py')),
        'verifier_sha256': sha(Path(__file__)),
    }
    record_path.write_text(json.dumps(result, indent=2) + '\n')
    record_path.with_suffix('.progress.json').write_text(json.dumps(
        {'status': 'complete', 'checked_rows': len(entries),
         'negative_rows': checked_negative,
         'elapsed_seconds': result['elapsed_seconds']}) + '\n')
    print(json.dumps(result))


def main():
    sys.set_int_max_str_digits(20000)
    parser = argparse.ArgumentParser()
    for name in ('balls', 'masks', 'certificate', 'galois', 'record'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    check(args.balls, args.masks, args.certificate, args.galois, args.record)


if __name__ == '__main__':
    main()
