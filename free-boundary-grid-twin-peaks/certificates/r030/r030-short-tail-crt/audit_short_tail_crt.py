#!/usr/bin/env python3
"""Independent CRT/Pratt audit of recovered n=7 short-tail residues."""
from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
from hashlib import sha256
from math import factorial, gcd, prod
from pathlib import Path
import json
import re

REQUIRED_SOURCE = '9d757fcfd65729242d2742cc94b9c34e3e88ba23d1052b985c3b70f0b8a02d13'
getcontext().prec = 60


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def verify_pratt(data: dict) -> list[int]:
    nodes = data['nodes']
    done: set[int] = set()
    def check(p: int) -> None:
        if p in done:
            return
        node = nodes.get(str(p))
        assert node is not None and int(node['p']) == p
        if p == 2:
            assert node == {'p': '2', 'base': 'two'}
            done.add(p)
            return
        factors = {int(q): int(e) for q, e in node['p_minus_1_factorization'].items()}
        assert p > 2 and factors
        value = 1
        for q, e in factors.items():
            assert q >= 2 and e >= 1
            check(q)
            value *= q ** e
        assert value == p - 1
        a = int(node['witness'])
        assert 1 < a < p and pow(a, p - 1, p) == 1
        for q in factors:
            assert gcd(pow(a, (p - 1) // q, p) - 1, p) == 1
        done.add(p)
    roots = [int(x) for x in data['roots']]
    for root in roots:
        check(root)
    return roots, len(done)


def crt_direct(residues: list[int], primes: list[int]) -> int:
    modulus = prod(primes)
    answer = 0
    for residue, prime in zip(residues, primes):
        cofactor = modulus // prime
        answer = (answer + residue * cofactor * pow(cofactor, -1, prime)) % modulus
    return answer


def crt_garner(residues: list[int], primes: list[int]) -> int:
    answer, modulus = 0, 1
    for residue, prime in zip(residues, primes):
        t = ((residue - answer) % prime) * pow(modulus, -1, prime) % prime
        answer += modulus * t
        modulus *= prime
    assert 0 <= answer < modulus
    return answer


def fdata(x: int, y: int) -> dict:
    q = Fraction(x, y)
    return {'numerator': str(q.numerator), 'denominator': str(q.denominator),
            'display': str(Decimal(q.numerator) / Decimal(q.denominator))}


def main(raw, pratt, reference, r29_source, output) -> None:
    cert = json.loads(pratt.read_text())
    primes, node_count = verify_pratt(cert)
    assert len(primes) == 4 and len(set(primes)) == 4
    records = []
    for pi, prime in enumerate(primes):
        json_path = raw / f'n7-v24-p{pi}.json'
        exit_path = raw / f'n7-v24-p{pi}.exit-status.txt'
        resource_path = raw / f'n7-v24-p{pi}.resource.txt'
        row = json.loads(json_path.read_text())
        assert exit_path.read_text().strip() == '0'
        resource = resource_path.read_text()
        assert re.search(r'\n\s*0\s+swaps\n', resource)
        assert row['n'] == 7 and row['v'] == 24 and row['prime_index'] == pi
        assert int(row['prime']) == prime and row['algorithm'] == 'proposed-short-tail-modular-layer-rolling'
        assert row['formula_status'] == 'pending independent proof'
        assert row['independent_sets'] == 293815906 and row['max_layer_count'] == 56669220
        records.append((row, json_path, exit_path, resource_path))
    raw_source = raw / 'short_tail_mod.cpp'
    assert digest(raw_source) == REQUIRED_SOURCE
    assert digest(r29_source) == REQUIRED_SOURCE
    modulus = prod(primes)
    bound = factorial(48)
    assert modulus > bound
    values = {}
    for key in ('h1', 'h2', 'h3'):
        residues = [int(row[f'{key}_residue']) for row, *_ in records]
        a = crt_direct(residues, primes)
        b = crt_garner(residues, primes)
        assert a == b
        assert 0 <= a <= bound
        assert all(a % p == r for p, r in zip(primes, residues))
        values[key] = a
    ref = json.loads(reference.read_text())
    assert ref['verified'] and ref['n'] == 7 and ref['v'] == 24
    A, B = int(ref['A2v']), int(ref['Bv'])
    short_mass = sum(values.values())
    short_moment = values['h1'] + 2 * values['h2'] + 3 * values['h3']
    remainder_mass = B - short_mass
    remainder_moment = A - short_moment
    assert remainder_mass >= 0 and remainder_moment >= 4 * remainder_mass and remainder_moment <= 48 * remainder_mass
    out = {
        'status': 'passed',
        'scope': 'exact finite n=7 center short-tail recovery; no asymptotic inference',
        'n': 7, 'v': 24,
        'source_sha256_verified': REQUIRED_SOURCE,
        'raw_formula_status_preserved': 'pending independent proof',
        'raw_formula_status_explanation': 'The Mac checkpoint labels are historical immutable metadata. This audit applies the subsequently accepted r29 identity but does not rewrite those raw labels.',
        'pratt': {'certificate_path': pratt.name, 'certificate_sha256': digest(pratt), 'roots': [str(x) for x in primes], 'verified_nodes': node_count},
        'crt': {'algorithms': ['direct cofactor CRT', 'incremental Garner CRT'], 'modulus': str(modulus), 'modulus_bits': modulus.bit_length(),
                'count_upper_bound': str(bound), 'count_upper_bound_bits': bound.bit_length(), 'unique_under_bound': modulus > bound},
        'h': {key: str(value) for key, value in values.items()},
        'reference_selected_normalizer': {'path': reference.name, 'sha256': digest(reference), 'A2v': str(A), 'Bv': str(B)},
        'short_mass': str(short_mass), 'short_moment': str(short_moment),
        'short_mass_over_Bv': fdata(short_mass, B), 'short_moment_over_A2v': fdata(short_moment, A),
        'residual_L_at_least_4_mass': str(remainder_mass),
        'residual_L_at_least_4_first_moment': str(remainder_moment),
        'residual_mass_over_Bv': fdata(remainder_mass, B),
        'residual_conditional_mean_L': fdata(remainder_moment, remainder_mass),
        'residual_constraints': {'support': '4 <= L <= 48', 'four_R_le_AminusS_le_48_R': True,
                                 'four_R': str(4 * remainder_mass), 'AminusS': str(remainder_moment), '48_R': str(48 * remainder_mass)},
        'raw_files': [{'json': j.name, 'json_sha256': digest(j), 'exit': e.name, 'exit_sha256': digest(e),
                       'resource': r.name, 'resource_sha256': digest(r)} for _, j, e, r in records]
    }
    output.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser(); p.add_argument('--raw',type=Path,required=True); p.add_argument('--pratt',type=Path,required=True); p.add_argument('--reference',type=Path,required=True); p.add_argument('--r29-source',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); main(a.raw,a.pratt,a.reference,a.r29_source,a.output)
