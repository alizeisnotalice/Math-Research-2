#!/usr/bin/env python3
"""Read-only replicate-weight diagnostic for completed n=24 terminal raws.

This is not a genealogy analysis: terminal records do not retain parent IDs,
offspring counts, or resampling ancestry.  It distinguishes the recorded
across-replicate normalizer Z from the recorded L-lift multiplier only.
"""
from __future__ import annotations

import hashlib
import json
import math
import pathlib
import statistics
import subprocess
import sys
from decimal import Decimal, getcontext

getcontext().prec = 90
BASE = 0xA022000000000000
ROOT = pathlib.Path(__file__).resolve().parent
ARCHIVE = ROOT.parent
FROZEN_AGG = ARCHIVE / 'r025-local-n24' / 'aggregate_tau2_distance.py'
FROZEN_SRC = ARCHIVE / 'r025-local-n24' / 'tau2_distance_terminal_r023.swift'

INPUTS = {
    'local_strat': ARCHIVE / 'r025-local-n24' / 'n24-p131072-r64-l16.jsonl',
    'macbook_multi': ARCHIVE / 'r025-macbook-gpu' / 'n24-distance-p131072-r64-l16.jsonl',
}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def d(x):
    return Decimal(str(x))

def pearson(xs, ys):
    if len(xs) < 2:
        return None
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    vx = sum((x-mx)**2 for x in xs)
    vy = sum((y-my)**2 for y in ys)
    return None if not vx or not vy else sum((x-mx)*(y-my) for x, y in zip(xs, ys)) / math.sqrt(vx*vy)

def normalize(ws):
    total = sum(ws)
    assert total > 0
    return [w/total for w in ws]

def ess(ps):
    return float(Decimal(1) / sum(p*p for p in ps))

def sidecar_for(raw, meta):
    expected = {k: meta[k] for k in ('n','particles','replicates','tau_lo','tau_hi','domain','method','move_steps')}
    p = ROOT / (raw.stem + '.sidecar.json')
    p.write_text(json.dumps({'raw_sha256': sha(raw), 'source_sha256': sha(FROZEN_SRC),
                             'expected_metadata': expected}, indent=2) + '\n')
    return p

def frozen_summary(raw, side):
    cmd = [sys.executable, str(FROZEN_AGG), str(raw), '--sidecar', str(side),
           '--source', str(FROZEN_SRC), '--epsilon', '1/2']
    run = subprocess.run(cmd, text=True, capture_output=True)
    if run.returncode:
        raise RuntimeError(run.stderr)
    return json.loads(run.stdout)

def analyze(label, raw):
    lines = [json.loads(x) for x in raw.read_text().splitlines() if x.strip()]
    assert len(lines) >= 3 and lines[0].get('kind') == 'metadata' and lines[-1].get('kind') == 'run_end'
    meta, records = lines[0], lines[1:-1]
    n, P, R = meta['n'], meta['particles'], meta['replicates']
    assert n == 24 and P == 131072 and R == 64 and meta['tau_lo'] == meta['tau_hi'] == 2
    assert len(records) == R and sorted(x['replicate'] for x in records) == list(range(R))
    assert all(x['kind'] == 'replicate' and x['tau'] == 2 and not x['dead'] for x in records)
    for rec in records:
        want = (BASE + int(meta['domain'], 16) + n*1_000_000 + 2*10_000 + rec['replicate']) & ((1 << 64)-1)
        assert rec['seed_hex'].lower().zfill(16) == f'{want:016x}'
    side = sidecar_for(raw, meta)
    frozen = frozen_summary(raw, side)

    logs = [d(x['logZ']) for x in records]
    maxlog = max(logs)
    z = [(x-maxlog).exp() for x in logs]
    # Reconstruct from integer group payloads, matching the frozen reducer
    # rather than trusting separately rounded producer summary fields.
    L = [sum(map(Decimal, x['tau2_lift_L_groups'])) / P for x in records]
    DD = [sum(map(Decimal, x['tau2_lift_L_distance2_groups'])) / P for x in records]
    hist = [[sum(Decimal(row[j]) for row in x['tau2_distance_L_hist_groups']) / P
             for j in range(2*n-1)] for x in records]
    assert all(len(v) == 2*n-1 for v in hist)
    zprob = normalize(z)
    lift = [a*b for a, b in zip(z, L)]
    lprob = normalize(lift)
    denom = sum(lift)
    q = sum(a*b for a, b in zip(z, DD)) / (Decimal(n*n)*denom)
    histfull = [sum(z[r]*hist[r][j] for r in range(R))/denom for j in range(2*n-1)]
    assert abs(sum(histfull)-1) < Decimal('1e-70')
    mean_manhattan = sum(Decimal(j)*x for j, x in enumerate(histfull))
    each_q = [DD[i]/(Decimal(n*n)*L[i]) for i in range(R)]
    each_mean_d = [sum(Decimal(j)*hist[i][j] for j in range(2*n-1))/L[i] for i in range(R)]
    loo = []
    for omit in range(R):
        d0 = denom-lift[omit]
        loo.append(sum(z[i]*DD[i] for i in range(R) if i != omit)/(Decimal(n*n)*d0))
    # These Pearson statistics are descriptive float calculations, not tests.
    flogs = [float(x) for x in logs]
    relation = {
        'corr_logZ_L_mean': pearson(flogs, [float(x) for x in L]),
        'corr_logZ_conditional_Q': pearson(flogs, [float(x) for x in each_q]),
        'corr_logZ_conditional_mean_Manhattan_distance': pearson(flogs, [float(x) for x in each_mean_d]),
        'corr_logZ_log_L': pearson(flogs, [math.log(float(x)) for x in L]),
        'corr_log_lift_weight_conditional_Q': pearson([math.log(float(a))+math.log(float(b)) for a,b in zip(z,L)], [float(x) for x in each_q]),
    }
    def weights(ps):
        i = max(range(R), key=lambda j: ps[j])
        return {'ESS': ess(ps), 'max_share': float(ps[i]), 'max_replicate': i,
                'top5_shares': [float(x) for x in sorted(ps, reverse=True)[:5]]}
    out = {
        'raw_sha256': sha(raw), 'input_validation': {
            'metadata': {k: meta[k] for k in ('n','particles','replicates','domain','method','move_steps')},
            'metadata_then_64_replicates_then_run_end': True,
            'deterministic_seeds_complete_unique': True,
            'all_tau2_records_nondead': True,
            'frozen_aggregator_accepted': True,
            'frozen_aggregator_sha256': sha(FROZEN_AGG), 'frozen_source_sha256': sha(FROZEN_SRC),
        },
        'frozen_aggregate_crosscheck': {
            'Q': frozen['Q'], 'tau2_ESS': frozen['tau2_weight']['ESS'],
            'lift_ESS': frozen['liftZ_weight']['ESS'], 'lift_max_share': frozen['liftZ_weight']['maxshare'],
        },
        'terminal_weight_diagnostic': {
            'logZ_min_mean_max_sd': [str(min(logs)), str(sum(logs)/R), str(max(logs)), str(statistics.stdev(map(float,logs)))],
            'E_tau2_L': str(sum(z[i]*L[i] for i in range(R))/sum(z)),
            'unweighted_L_min_mean_max': [str(min(L)), str(sum(L)/R), str(max(L))],
            'Q_complete_K2_L_lift_ratio': str(q),
            'target_mean_Manhattan_distance': str(mean_manhattan),
            'complete_K2_L_weighted_distance_histogram': [str(x) for x in histfull],
            'tau2_normalizer_weights': weights(zprob),
            'complete_K2_L_lift_weights': weights(lprob),
            'leave_one_replicate_out_Q_min_max': [str(min(loo)), str(max(loo))],
            'leave_one_replicate_out_Q_at_full_weight_max': str(loo[max(range(R), key=lambda j:lprob[j])]),
            'per_replicate_conditional_Q_min_mean_max': [str(min(each_q)),str(sum(each_q)/R),str(max(each_q))],
            'logZ_geometry_descriptive_correlations': relation,
        },
        'interpretation_scope': 'Terminal replicate-normalizer diagnostic only. No parent IDs, offspring counts, ESS-by-rank, resampling times, terminal ancestry, or trajectory genealogy are retained, so no genealogy claim can be made.',
        'no_pooling': 'local strat and MacBook multi use different proposals and random domains; reported separately.',
        'uncertainty_scope': 'No delta-SE is interpreted as a confidence interval or mixing certificate.'
    }
    assert abs(q-d(frozen['Q'])) < Decimal('1e-70')
    return out

out = {'kind': 'n24_terminal_weight_diagnostic', 'runs': {label: analyze(label, raw) for label, raw in INPUTS.items()}}
(ROOT/'n24-genealogy-diagnostic.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
