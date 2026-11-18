#!/usr/bin/env python3
"""Independently transform the real restricted truth tables in cut obstructions."""
from collections import Counter
from fractions import Fraction
import json
import platform
import time

import scan_cut as scan


def witnesses():
    return [
        ('six-bit-explicit-r49',6,scan.fixture_truth(6,lambda u:scan.r49((u&7)|((u&~7)<<1)))),
        ('six-bit-macbook-cubic',6,1599346001363546751),
        ('six-bit-mini-cubic',6,11679430202935411900),
        ('seven-bit-r49',7,scan.fixture_truth(7,scan.r49)),
        ('seven-bit-weak-folding',7,scan.fixture_truth(7,scan.weak_fold_1)),
    ]


def direct_child(truth,n,direction,outcome):
    lead=1<<(direction.bit_length()-1)
    reps=scan.representatives(n,lead)
    values=[]
    for z in reps:
        physical=z|((scan.parity(z&direction)^outcome)*lead)
        values.append(1-2*((truth>>physical)&1))
    return {u:w for u,w in zip(reps,scan.fwht(values)) if w}


def anf_terms(truth,n):
    bits=[(truth>>u)&1 for u in range(1<<n)]
    for i in range(n):
        for u in range(1<<n):
            if u&(1<<i):
                bits[u]^=bits[u^(1<<i)]
    return [u for u,b in enumerate(bits) if b]


def verify(name,n,truth):
    record=scan.analyze(truth,n,True)
    histogram=Counter()
    max_child_sq=Fraction(0)
    for row in record['all_cuts']:
        q=row['query']
        left=direct_child(truth,n,q,0)
        right=direct_child(truth,n,q,1)
        if left!=row['walsh0'] or right!=row['walsh1']:
            raise AssertionError((name,q,'direct restricted truth disagreement'))
        if not row['effective']:
            continue
        if row['signed_margin']>=0 or row['sharpened_margin']>=0 or max(row['child_vs_parent_numerators'])>=0:
            raise AssertionError((name,q,'not an all-effective-directions obstruction'))
        max_child_sq=max(max_child_sq,*(Fraction((row[f'R{i}']+1)**2,row[f's{i}']) for i in (0,1)))
        key=tuple(row[k] for k in ('s0','s1','R0','R1','j','k','e','signed_margin'))
        histogram[key]+=1
    return dict(name=name,anf_terms=anf_terms(truth,n),
                best_child_P_squared=[max_child_sq.numerator,max_child_sq.denominator],
                checked_restricted_truth_tables=2*((1<<n)-1),
                cut_class_histogram=[dict(zip(('s0','s1','R0','R1','j','k','e','signed_margin','directions'),(*k,v)))
                                     for k,v in sorted(histogram.items())],**record)


if __name__=='__main__':
    start=time.perf_counter()
    records=[verify(*w) for w in witnesses()]
    print(json.dumps(dict(status='pass',wall_seconds=time.perf_counter()-start,host=platform.node(),
                          python=platform.python_version(),records=records),indent=2,sort_keys=True))
