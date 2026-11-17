#!/usr/bin/env python3
"""Independent reconstruction of two selected-hole counts, not full A2(7)."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('certificate', type=Path)
    ap.add_argument('residues', type=Path, nargs=4)
    args = ap.parse_args()
    cert = json.loads(args.certificate.read_text())
    nodes, done = cert['nodes'], set()
    def prime(p):
        if p in done:
            return
        node = nodes[str(p)]
        assert int(node['p']) == p
        if p == 2:
            assert node == {'p':'2', 'base':'two'}
            done.add(p)
            return
        factors = {int(q):int(e) for q,e in node['p_minus_1_factorization'].items()}
        assert factors and all(2 <= q < p and e >= 1 for q,e in factors.items())
        assert math.prod(q**e for q,e in factors.items()) == p-1
        for q in factors:
            prime(q)
        a = int(node['witness'])
        assert 1 < a < p and pow(a,p-1,p) == 1
        assert all(math.gcd(pow(a,(p-1)//q,p)-1,p) == 1 for q in factors)
        done.add(p)
    primes = list(map(int,cert['roots']))
    assert len(primes) == len(set(primes)) == 4
    for p in primes:
        prime(p)
    # Independent forward row transfer with a hole and its four neighbors removed.
    n, v = 7, 24
    forbidden = {v,v-n,v+n,v-1,v+1}
    dp = {(0,0):1}
    for y in range(n):
        blocked = sum(1 << x for x in range(n) if y*n+x in forbidden)
        valid = [r for r in range(1<<n) if not r&(r<<1) and not r&blocked]
        nxt = {}
        for (prev,k),count in dp.items():
            for r in valid:
                if not prev&r:
                    key = (r,k+r.bit_count())
                    nxt[key] = nxt.get(key,0)+count
        dp = nxt
    histogram = [0]* (1+max(k for _,k in dp))
    for (_,k),count in dp.items():
        histogram[k] += count
    assert sum(histogram) == 293815906 and max(histogram) == 56669220
    rows = sorted([json.loads(p.read_text()) for p in args.residues],key=lambda r:r['prime_index'])
    assert [r['prime_index'] for r in rows] == list(range(4))
    for i,r in enumerate(rows):
        assert r['n'] == n and r['v'] == v
        assert r['algorithm'] == 'punctured-modular-layer-rolling'
        assert int(r['prime']) == primes[i] > 48
        assert r['independent_sets'] == sum(histogram)
        assert r['max_layer_count'] == max(histogram)
        assert r['max_modeled_payload_bytes'] <= r['payload_limit_bytes']
    upper = math.factorial(48)
    M = math.prod(primes)
    assert M > upper
    values = {}
    for key in ('A2v_residue','Bv_residue'):
        residues = [int(r[key]) for r in rows]
        assert all(0 <= a < p for a,p in zip(residues,primes))
        x,m = 0,1
        for a,p in zip(residues,primes):
            x += m*((a-x)*pow(m,-1,p)%p)
            m *= p
        y = sum(a*(M//p)*pow(M//p,-1,p) for a,p in zip(residues,primes)) % M
        assert x == y and all(x % p == a for a,p in zip(residues,primes))
        assert 0 < x <= upper
        values[key.replace('_residue','')] = x
    A,B = values['A2v'],values['Bv']
    assert B <= A
    delta = 5*B-A
    out = {'verified':True,'n':n,'v':v,'A2v':str(A),'Bv':str(B),
        'scope':'exact selected central-hole counts; no asymptotic claim',
        'Bv_over_A2v':{'numerator':str(B),'denominator':str(A),'display':B/A},
        'five_Bv_minus_A2v':str(delta),'Bv_over_A2v_at_least_one_fifth':delta>=0,
        'CRT_modulus':str(M),'upper_bound_48_factorial':str(upper),
        'CRT_unique_within_proved_nonnegative_count_bound':True,
        'primality_nodes_verified':len(done),'independent_set_layer_histogram':histogram,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.certificate,*args.residues]}}
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
