"""Generate recursive Lucas/Pratt certificates for four 61-bit candidates."""
from collections import Counter
from math import gcd, prod
from pathlib import Path
import json


def trial(n):
    """A candidate filter only. The proof is the recursive certificate."""
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    odd = n - 1
    shifts = 0
    while odd % 2 == 0:
        odd //= 2
        shifts += 1
    for b in (2, 3, 5, 7, 11, 13, 17, 19, 23):
        x = pow(b, odd, n)
        if x in (1, n - 1):
            continue
        for _ in range(shifts - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factor(n):
    if n == 1:
        return []
    if trial(n):
        return [n]
    for q in range(2, 100):
        if n % q == 0:
            return factor(q) + factor(n // q)
    for c in range(1, 100):
        x = y = 2
        d = 1
        while d == 1:
            x = (x*x + c) % n
            y = (y*y + c) % n
            y = (y*y + c) % n
            d = gcd(abs(x-y), n)
        if d < n:
            return factor(d) + factor(n // d)
    raise ValueError(f'factorization failed for {n}')


def build(p, tree):
    if str(p) in tree:
        return
    if p == 2:
        tree['2'] = {'p': '2', 'base': 'two'}
        return
    factors = Counter(factor(p-1))
    assert prod(q**e for q, e in factors.items()) == p-1
    for q in factors:
        build(q, tree)
    a = next(a for a in range(2, 5000) if pow(a, p-1, p) == 1
             and all(gcd(pow(a, (p-1)//q, p)-1, p) == 1 for q in factors))
    tree[str(p)] = {'p': str(p), 'witness': str(a),
                    'p_minus_1_factorization': {str(q): e for q, e in sorted(factors.items())}}


def main():
    choices = []
    proof = {}
    for p in range(2**61-1, 2**61-10000, -2):
        if trial(p):
            build(p, proof)
            choices.append(str(p))
        if len(choices) == 4:
            break
    assert len(choices) == 4
    dest = Path(__file__).with_name('pratt-61bit.json')
    dest.write_text(json.dumps({'roots': choices, 'nodes': proof}, indent=2, sort_keys=True) + '\n')
    print(', '.join(choices))


if __name__ == '__main__':
    main()
