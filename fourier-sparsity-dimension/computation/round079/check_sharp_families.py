#!/usr/bin/env python3
"""Verify six-point Boolean fillings and sharp multi-circuit Fourier flats."""
from collections import Counter
import hashlib
import json
import platform
import time

from check_cut import gf2_basis, parity, six_circuit, walsh


def fill_census():
    fixture = six_circuit()
    g = fixture['ternary_values']
    zero = [x for x,u in enumerate(g) if not u]
    censored = Counter()
    witnesses = {}
    for mask in range(1 << len(zero)):
        h = {x:1-2*((mask>>i)&1) for i,x in enumerate(zero)}
        values = [a or (1-2*z)*h[x] for z in (0,1) for x,a in enumerate(g)]
        truth = sum((a == -1) << x for x,a in enumerate(values))
        raw = walsh(truth,5)
        if raw[:16] != fixture['full_integer_walsh'][:16]:
            raise AssertionError('fill changed the ternary flat')
        support = [u for u,a in enumerate(raw) if a]
        d = len(gf2_basis(support))
        censored[len(support),d] += 1
        witnesses.setdefault((len(support),d), dict(truth_mask=truth,full_integer_walsh=raw))
    return dict(fill_count=1 << len(zero),zero_points=zero,
                spectra=[dict(s=s,d=d,count=count,example=witnesses[s,d])
                         for (s,d),count in sorted(censored.items())])


def sharp_case(t, lengths):
    ell = len(lengths)
    if any(m < 4 or m % 2 for m in lengths) or sum(m-2 for m in lengths) != 1 << (t + 1):
        raise ValueError('circuit lengths give the wrong extreme lattice range')
    spectrum = {}
    next_offset = ell if ell > 1 else 0
    for block,m in enumerate(lengths):
        tag = 1 << block if ell > 1 else 0
        offsets = [1 << (next_offset + k) for k in range(m - 2)]
        next_offset += m - 2
        for point in [0,*offsets,sum(offsets)]:
            spectrum[tag ^ point] = -1 if point == 0 else 1
    n = next_offset
    count = 1 << (t+n)
    denom = 1 << (t+1)
    lattice = Counter()
    truth_mask = 0
    for x in range(1 << n):
        numerator = sum(a*(1-2*parity(x&u)) for u,a in spectrum.items())
        if not -denom <= numerator <= denom:
            raise AssertionError('lattice value out of range')
        negatives,residue = divmod(denom - numerator, 4)
        if residue or not 0 <= negatives <= 1 << t:
            raise AssertionError('lattice value not realizable by bit filling')
        lattice[numerator] += 1
        for y in range(negatives):
            truth_mask |= 1 << (x + (y << n))
    full = walsh(truth_mask,n+t)
    base = full[:(1 << n)]
    if any(value*denom != spectrum.get(u,0)*count for u,value in enumerate(base)):
        raise AssertionError('full Boolean transform does not contain the claimed flat')
    q = sum(a != 0 for a in base)
    rank = len(gf2_basis((1 << n)|u for u,a in enumerate(base) if a))
    if q-rank != ell or q != 2*ell + (1 << (t+1)):
        raise AssertionError('actual Boolean spectrum is not an equality case')
    return dict(t=t,ell=ell,ambient_dimension=n+t,lengths=lengths,q=q,affine_rank=rank,
                full_sparsity=sum(a!=0 for a in full),lattice_numerator_histogram=dict(lattice),
                canonical_truth_sha256=hashlib.sha256(truth_mask.to_bytes((count+7)//8,'little')).hexdigest())


def main():
    start, cpu = time.perf_counter(),time.process_time()
    cases = [(0,[4]),(1,[6]),(1,[4,4]),(2,[10]),(2,[8,4]),(2,[6,4,4]),(2,[4,4,4,4]),(3,[18])]
    records = [sharp_case(t,m) for t,m in cases]
    print(json.dumps(dict(status='pass',python=platform.python_version(),host=platform.node(),
                          six_circuit_fillings=fill_census(),sharp_examples=records,
                          wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time()-cpu),
                     indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
