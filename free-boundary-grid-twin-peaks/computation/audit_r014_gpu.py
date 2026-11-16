#!/usr/bin/env python3
"""Audit finite n8/n9 SMC data, with Decimal arithmetic and explicit pooling.

All reported standard errors are empirical delta diagnostics. No particle is
treated as an independent repeat. PRNG seeds and finite repetitions do not
certify independence, bias, confidence coverage, or an asymptotic statement.
"""
from collections import defaultdict
from decimal import Decimal as D, getcontext
from pathlib import Path
import csv
import hashlib
import json
import math
from audit_tau_smc_aggregates import ROOT, aggregate, read

getcontext().prec = 60
BASE = ROOT / 'certificates/r014-gpu'


def meanvar(values):
    n = len(values)
    assert n > 1
    mean = sum(values) / n
    return mean, sum((x-mean)**2 for x in values) / (n*(n-1))


def metal(path):
    return path.read_text().split('let metal = """', 1)[1].split('"""', 1)[0]


def audit():
    seen = set()
    source = ROOT / 'certificates/r013-gpu/audit/tau_smc_r013.swift'
    assert metal(source) == metal(BASE/'production_r013.swift')
    assert metal(source) == metal(BASE/'onepeak/source/onepeak_smc.swift')
    two_rows, two_data, telemetry = [], {}, []
    for path in sorted(BASE.glob('n*.jsonl')):
        m, groups, raw = read(path)
        assert raw[-1]['kind'] == 'run_end'
        assert m['tau_lo'] == 2 and m['tau_hi'] == m['n']**2
        for t, reps in groups.items():
            for r, x in reps.items():
                seed = (0xA022000000000000 + int(m['seed_domain_hex'],16)
                        + m['n']*1000000 + t*10000 + r) % (1 << 64)
                assert int(x['seed_hex'], 16) == seed and seed not in seen
                assert x['dead'] or x['Z'] > 0
                seen.add(seed)
        row, _ = aggregate(path)
        two_rows.append(row)
        two_data[path.name] = (m, groups)
        monitor = path.with_name(path.stem+'-monitor.csv')
        if monitor.exists():
            records = list(csv.DictReader(monitor.open()))
            vals = [D(x['agx_device_util_pct']) for x in records if x.get('agx_device_util_pct')]
            telemetry.append({'path':str(monitor.relative_to(ROOT)),
                              'samples_with_device_reading':len(vals),
                              'device_mean_percent':str(sum(vals)/len(vals)) if vals else None,
                              'samples_at_100':sum(x == 100 for x in vals),
                              'nonexclusive_whole_device':True})
    one_rows, one_data = [], defaultdict(list)
    for path in sorted((BASE/'onepeak/raw').glob('*.jsonl')):
        raw = [json.loads(s,parse_float=D) for s in path.read_text().splitlines()]
        m = raw[0]
        assert m['kind'] == 'metadata' and raw[-1]['kind'] == 'run_end'
        assert m['tau'] == m['n']**2 + 1
        reps = [x for x in raw if x['kind'] == 'replicate']
        assert [x['replicate'] for x in reps] == list(range(m['replicates']))
        for r, x in enumerate(reps):
            seed = (int(m['seed_domain_hex'],16)+m['n']*1000000+r) % (1 << 64)
            assert int(x['seed_hex'],16) == seed and seed not in seen
            assert x['tau'] == m['tau'] and not x['dead'] and 0 < x['Z'] <= 1
            seen.add(seed)
        z, variance = meanvar([x['Z'] for x in reps])
        one_rows.append({'path':str(path.relative_to(ROOT)), 'n':m['n'],
                         'particles':m['particles'],'repeats':len(reps),
                         'Z':str(z),'empirical_SE':str(variance.sqrt()),
                         'dead':0,'wall_seconds':str(raw[-1]['wall_seconds'])})
        if m['n'] in (8,9):
            assert m['particles'] == 131072
            one_data[m['n']].extend(x['Z'] for x in reps)
        if m['n'] == 7:
            exact = D(134267527713152677457983522684235911142268448963584)/D(math.factorial(49))
            one_rows[-1]['exact_Z'] = str(exact)
            one_rows[-1]['relative_error_to_exact'] = str(z/exact-1)
    # Pool only the same particle count/method. The r013 32-repeat n8 batch
    # is an explicitly reused dependency, not a new r014 observation.
    n8_old = ROOT/'certificates/r013-gpu/stratified/n8-p131072-r32-strat.jsonl'
    old_m, old_g, old_raw = read(n8_old)
    assert old_raw[-1]['kind'] == 'run_end'
    for t,g in old_g.items():
        for r,x in g.items():
            seed = (0xA022000000000000+0x7600000000000000+8*1000000+t*10000+r) % (1<<64)
            assert int(x['seed_hex'],16)==seed and seed not in seen
            seen.add(seed)
    pooled = []
    for n in (8,9):
        named = [(name,m,g) for name,(m,g) in two_data.items()
                 if m['n'] == n and m['particles'] == 131072]
        if n == 8:
            named.append((str(n8_old.relative_to(ROOT)),old_m,old_g))
        strata = defaultdict(list)
        for _,m,groups in named:
            for t,g in groups.items(): strata[t].extend(g.values())
        z2 = sum(meanvar([x['Z'] for x in g])[0] for g in strata.values())
        u2 = sum(meanvar([x['U'] for x in g])[0] for g in strata.values())
        v2 = sum(meanvar([x['Z'] for x in g])[1] for g in strata.values())
        q = u2/z2
        vq = sum(meanvar([x['U']-q*x['Z'] for x in g])[1] for g in strata.values())
        z1,v1 = meanvar(one_data[n])
        ratio = z2/z1
        vr = v2/z1**2 + z2**2*v1/z1**4
        pooled.append({'n':n,'particles':131072,'twopeak_inputs':[name for name,_,_ in named],
                       'twopeak_repeats_per_tau':len(strata[2]),'onepeak_repeats':len(one_data[n]),
                       'Z1':str(z1),'Z1_empirical_SE':str(v1.sqrt()),
                       'Z2':str(z2),'Z2_empirical_SE':str(v2.sqrt()),
                       'Q':str(q),'Q_empirical_delta_SE':str(vq.sqrt()/z2),
                       'R':str(ratio),'R_empirical_delta_SE':str(vr.sqrt()),
                       'R_over_N':str(ratio/D(n*n)),
                       'logR_over_nlogn':str(ratio.ln()/(D(n)*D(n).ln()))})
    sources = [{'path':str(p.relative_to(ROOT)), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
               for p in (source,BASE/'production_r013.swift',BASE/'onepeak/source/onepeak_smc.swift')]
    return {'status':'passed','precision_decimal_digits':60,'new_and_reused_seeds_unique':len(seen),
            'sources':sources,'twopeak_batches':two_rows,'onepeak_batches':one_rows,
            'same_P_pooled_diagnostics':pooled,'telemetry':telemetry,
            'uncertainties_are_descriptive_not_certificates':True}


if __name__ == '__main__':
    print(json.dumps(audit(),indent=2))
