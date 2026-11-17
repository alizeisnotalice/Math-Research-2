#!/usr/bin/env python3
"""Exact Walsh verification of bentness for the CDDL Fact-22 seed.

This checks only the bent assertion.  It does not, and cannot by itself,
replace CDDL's normality-search certificate for non-weak-normality.
"""

from math import gcd

N = 14
MASK = (1 << N) - 1


def degree(a: int) -> int:
    return a.bit_length() - 1


def poly_rem(a: int, modulus: int) -> int:
    d = degree(modulus)
    while a and degree(a) >= d:
        a ^= modulus << (degree(a) - d)
    return a


def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_rem(a, b)
    return a


def pmul(a: int, b: int, modulus: int) -> int:
    ans = 0
    while b:
        if b & 1:
            ans ^= a
        b >>= 1
        a <<= 1
        if a & (1 << N):
            a ^= modulus
    return ans & MASK


def ppow2(a: int, e: int, modulus: int) -> int:
    for _ in range(e):
        a = pmul(a, a, modulus)
    return a


def is_irreducible(modulus: int) -> bool:
    x = 2
    if ppow2(x, N, modulus) != x:
        return False
    # Prime divisors of 14 are 2 and 7.
    return all(poly_gcd(ppow2(x, N // q, modulus) ^ x, modulus) == 1 for q in (2, 7))


MOD = next(p for p in range((1 << N) | 1, 1 << (N + 1), 2) if is_irreducible(p))


def mul(a: int, b: int) -> int:
    return pmul(a, b, MOD)


def power(a: int, e: int) -> int:
    out = 1
    while e:
        if e & 1:
            out = mul(out, a)
        a = mul(a, a)
        e >>= 1
    return out


def trace(a: int) -> int:
    out = 0
    for _ in range(N):
        out ^= a
        a = mul(a, a)
    assert out in (0, 1)
    return out


def fwht(values: list[int]) -> None:
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for j in range(start, start + width):
                values[j], values[j + width] = values[j] + values[j + width], values[j] - values[j + width]
        width <<= 1


def main():
    import argparse,hashlib,json,sys,time
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();start=time.monotonic()
    alpha=next(x for x in range(2,1<<N) if mul(x,x)^x==1)
    assert power(alpha,3)==1
    truth=[trace(mul(alpha,power(x,57))) for x in range(1<<N)]
    walsh=[1-2*b for b in truth];fwht(walsh)
    hist={str(v):walsh.count(v) for v in sorted(set(walsh))}
    assert set(map(abs,walsh))=={128}
    result={'status':'passed','scope':'independent exact bent-only check; nonweak-normality relies on cited CDDL Fact22',
        'modulus_hex':hex(MOD),'alpha_hex':hex(alpha),'truth_weight':sum(truth),'walsh_histogram':hist,
        'truthmask_hex':hex(sum(b<<i for i,b in enumerate(truth))),
        'truth_bytes_sha256':hashlib.sha256(bytes(truth)).hexdigest(),
        'python':sys.version,'elapsed_seconds':time.monotonic()-start,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='truthmask_hex'}))
if __name__=='__main__':main()
