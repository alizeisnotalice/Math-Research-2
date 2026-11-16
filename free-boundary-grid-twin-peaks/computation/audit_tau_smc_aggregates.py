#!/usr/bin/env python3
"""Decimal reaggregation of r012 tau-SMC, including zero-weight repetitions.

Run from any cwd. The printed uncertainties are empirical delta diagnostics,
not confidence certificates and not variances of independent particles.
"""
from collections import defaultdict
from decimal import Decimal, getcontext
from pathlib import Path
import json

getcontext().prec = 60
ROOT = Path(__file__).resolve().parents[1]
D = Decimal


def read(path):
    rows = [json.loads(s, parse_float=D) for s in path.read_text().splitlines()]
    meta = next(x for x in rows if x.get('kind') == 'metadata')
    groups = defaultdict(dict)
    for x in rows:
        if x.get('kind') != 'replicate':
            continue
        t, r = x['tau'], x['replicate']
        assert r not in groups[t]
        assert 0 <= x['Z'] <= 1 and 0 <= x['U'] <= 2*x['Z']
        if x['dead']:
            assert x['Z'] == x['U'] == 0 and x['logZ'] is None
        groups[t][r] = x
    assert set(groups) == set(range(meta['tau_lo'], meta['tau_hi']+1))
    for g in groups.values():
        assert set(g) == set(range(meta['replicates']))
    return meta, groups, rows


def aggregate(path):
    m, g, rows = read(path)
    R = m['replicates']
    z = sum(sum(x['Z'] for x in a.values()) / R for a in g.values())
    u = sum(sum(x['U'] for x in a.values()) / R for a in g.values())
    q = u/z
    variance = D(0)
    for a in g.values():
        values = [x['U']-q*x['Z'] for x in a.values()]
        mean = sum(values)/R
        variance += sum((v-mean)**2 for v in values)/(R*(R-1))
    zs = [sum(a[r]['Z'] for a in g.values()) for r in range(R)]
    us = [sum(a[r]['U'] for a in g.values()) for r in range(R)]
    residuals = [b-q*a for a,b in zip(zs,us)]
    mean = sum(residuals)/R
    paired_var = sum((v-mean)**2 for v in residuals)/(R*(R-1))
    out = dict(path=str(path.relative_to(ROOT)), n=m['n'], particles=m['particles'],
               repetitions_per_tau=R, strata=len(g), Z=str(z), U=str(u), Q=str(q),
               stratified_delta_se=str(variance.sqrt()/z),
               paired_replicate_delta_se=str(paired_var.sqrt()/z),
               total_normalizer_replicate_ESS=str(sum(zs)**2/sum(v*v for v in zs)),
               zero_repetitions=sum(x['dead'] for a in g.values() for x in a.values()))
    for key in ('gpu_hardware_seconds', 'synchronous_dispatch_seconds'):
        if all(key in x for a in g.values() for x in a.values()):
            out[key] = str(sum(x[key] for a in g.values() for x in a.values()))
    end = [x for x in rows if x.get('kind') == 'run_end']
    if end:
        assert len(end) == 1
        out['wall_seconds'] = str(end[0]['wall_seconds'])
    return out, g


if __name__ == '__main__':
    output = []
    for n, name, newname in (
        (7, 'n7-p16384-r8-normalized.jsonl', 'n7-p16384-r8-new.jsonl'),
        (8, 'n8-p32768-r16-timed.jsonl', 'n8-p32768-r16-new.jsonl'),
    ):
        a, ga = aggregate(ROOT/'certificates/r012-tau-smc/data'/name)
        b, gb = aggregate(ROOT/'certificates/r012-gpu-parent/raw'/newname)
        for t in set(ga) | set(gb):
            if t not in ga or t not in gb:
                assert t == n*n
                assert all(x['Z'] == x['U'] == 0 for x in (ga.get(t) or gb[t]).values())
                continue
            for r in ga[t]:
                for key in ('logZ', 'Z', 'U', 'dead'):
                    assert ga[t][r][key] == gb[t][r][key], (n,t,r,key)
        assert a['Q'] == b['Q'] and a['Z'] == b['Z']
        output.extend([a,b])
    print(json.dumps({'status':'passed', 'same_seed_reruns_are_not_new_samples':True,
                      'rows':output}, indent=2))
