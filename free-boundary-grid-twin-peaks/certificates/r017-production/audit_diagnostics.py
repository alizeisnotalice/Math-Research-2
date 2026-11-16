#!/usr/bin/env python3
"""Read-only Decimal(80) audit of completed A022 production JSONL files."""
from collections import defaultdict
from decimal import Decimal, getcontext
import hashlib
import json
from pathlib import Path

getcontext().prec = 80
ROOT = Path(__file__).resolve().parent
MASK = (1 << 64) - 1
BASE = 0xA022000000000000

FILES = {
    'n16_all_p65536_d72': ROOT / 'n16-p65536-r16-domain72.jsonl',
    'n16_onepeak_p65536_d85': ROOT / 'n16-onepeak-p65536-r256-domain85.jsonl',
    'n16_onepeak_p65536_d86': ROOT / 'n16-onepeak-p65536-r256-domain86.jsonl',
    'n16_tau2_p262144_d73': ROOT / 'n16-tau2-p262144-r64-domain73.jsonl',
    'n12_p524288_d70': ROOT / 'n12-p524288-r16-domain70.jsonl',
    'n12_p262144_d68': ROOT.parent / 'r016-gpu/production/n12-p262144-r16-strat-domain6800.jsonl',
    'n12_p262144_d69': ROOT.parent / 'r016-gpu/production/n12-p262144-r16-strat-domain6900.jsonl',
}


def D(value):
    return Decimal(str(value))


def fmt(x):
    if x is None:
        return None
    return format(x, '.20E')


def affine_seed(domain, n, tau, replicate):
    return (BASE + domain + n * 1_000_000 + tau * 10_000 + replicate) & MASK


def sample_moments(xs, ys=None):
    m = len(xs)
    mean_x = sum(xs) / m
    variance_x = (sum((x - mean_x) ** 2 for x in xs) / (m - 1)) if m > 1 else Decimal(0)
    out = {'count': m, 'mean_Z': mean_x, 'sample_var_Z': variance_x,
           'var_mean_Z': variance_x / m}
    total = sum(xs)
    out['ESS_Z'] = (total * total / sum(x * x for x in xs)) if total else Decimal(0)
    out['maxshare_Z'] = max(xs) / total if total else Decimal(0)
    if ys is not None:
        mean_y = sum(ys) / m
        variance_y = (sum((y - mean_y) ** 2 for y in ys) / (m - 1)) if m > 1 else Decimal(0)
        cov = (sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / (m - 1)) if m > 1 else Decimal(0)
        out.update(mean_U=mean_y, sample_var_U=variance_y, sample_cov_ZU=cov,
                   var_mean_U=variance_y / m, cov_mean_ZU=cov / m)
    return out


def load(name, family):
    path = FILES[name]
    # Preserve each decimal JSON numeral directly; do not round through binary
    # float before the Decimal(80) aggregation.
    rows = [json.loads(line, parse_float=Decimal) for line in path.read_text().splitlines()]
    meta = rows[0]
    assert meta['kind'] == 'metadata' and rows[-1]['kind'] == 'run_end'
    n, p, r = meta['n'], meta['particles'], meta['replicates']
    domain_text = meta.get('domain', meta.get('seed_domain_hex'))
    domain = int(domain_text, 16)
    reps = [row for row in rows[1:-1] if row['kind'] == 'replicate']
    expected_count = (meta['tau_hi'] - meta['tau_lo'] + 1) * r
    assert len(reps) == expected_count, (name, len(reps), expected_count)
    by_tau = defaultdict(list)
    dead = 0
    for row in reps:
        expected = affine_seed(domain, n, row['tau'], row['replicate'])
        assert int(row['seed_hex'], 16) == expected, (name, row['tau'], row['replicate'])
        by_tau[row['tau']].append(row)
        if row['dead']:
            dead += 1
            assert D(row['Z']) == 0
            if row.get('U') is not None:
                assert D(row['U']) == 0
        else:
            # This particular completed n=12/n=16 input set has no live-Z
            # exponent underflow.  A future underflow-capable reader must use
            # logZ rather than apply this assertion universally.
            assert D(row['Z']) > 0 and row.get('logZ') is not None
    assert sorted(by_tau) == list(range(meta['tau_lo'], meta['tau_hi'] + 1))
    for tau, values in by_tau.items():
        assert sorted(row['replicate'] for row in values) == list(range(r)), (name, tau)
    return {'name': name, 'family': family, 'path': str(path.relative_to(ROOT.parent)),
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'meta': meta,
            'domain': domain, 'reps': reps, 'by_tau': by_tau, 'dead_count': dead,
            'coverage': {'expected_replicate_records': expected_count,
                         'actual_replicate_records': len(reps),
                         'all_tau_and_replicate_indices_complete': True,
                         'all_seed_formula_checks_passed': True}}


def per_tau_diagnostics(data):
    result = []
    for tau in sorted(data['by_tau']):
        rows = data['by_tau'][tau]
        zs = [D(row['Z']) for row in rows]  # dead rows deliberately remain zeros
        us = None if any(row.get('U') is None for row in rows) else [D(row['U']) for row in rows]
        moment = sample_moments(zs, us)
        result.append({
            'tau': tau, 'dead': sum(bool(row['dead']) for row in rows),
            'mean_Z': fmt(moment['mean_Z']), 'ESS_Z': fmt(moment['ESS_Z']),
            'maxshare_Z': fmt(moment['maxshare_Z']),
            'mean_U': fmt(moment.get('mean_U')),
        })
    return result


def summed_ratio(data, label):
    """Sum tau-specific means, with a stratified delta-method diagnostic SE."""
    zsum = usum = Decimal(0)
    var_z = var_u = cov_zu = Decimal(0)
    diag = []
    for tau in sorted(data['by_tau']):
        rows = data['by_tau'][tau]
        zs = [D(row['Z']) for row in rows]
        us = [D(row['U']) for row in rows]
        m = sample_moments(zs, us)
        zsum += m['mean_Z']; usum += m['mean_U']
        var_z += m['var_mean_Z']; var_u += m['var_mean_U']; cov_zu += m['cov_mean_ZU']
        diag.append({'tau': tau, 'ESS_Z': fmt(m['ESS_Z']), 'maxshare_Z': fmt(m['maxshare_Z']),
                     'dead': sum(bool(row['dead']) for row in rows)})
    q = usum / zsum
    var_q = var_u / (zsum * zsum) + usum * usum * var_z / (zsum ** 4) - Decimal(2) * usum * cov_zu / (zsum ** 3)
    # Decimal roundoff can make a mathematically nonnegative expression tiny negative.
    se = max(var_q, Decimal(0)).sqrt()
    T = len(data['by_tau'])
    uniform_z = zsum / T
    uniform_u = usum / T
    # The identity is algebraic; Decimal divisions round independently.
    assert abs(uniform_u / uniform_z - q) < Decimal('1e-70')
    return {
        'label': label, 'Z_sum_tau_means': fmt(zsum), 'U_sum_tau_means': fmt(usum),
        'Q_ratio': fmt(q), 'DeltaSE_diagnostic': fmt(se),
        'uniform_tau_mean_Z': fmt(uniform_z), 'uniform_tau_mean_U': fmt(uniform_u),
        'sum_vs_uniform_factor': str(T),
        'Q_uniform_ratio_identical': fmt(uniform_u / uniform_z),
        'per_tau_ESS_and_maxshare': diag,
        'warning': 'Diagnostic only: finite-R stratified delta calculation; it is not a certification or a coverage guarantee.'
    }


def combine_same_p(datas, label):
    """Pool only independent files with identical n/P/tau span, preserving tau strata."""
    first = datas[0]
    n, p = first['meta']['n'], first['meta']['particles']
    assert all(d['meta']['n'] == n and d['meta']['particles'] == p and sorted(d['by_tau']) == sorted(first['by_tau']) for d in datas)
    synthetic = {'by_tau': defaultdict(list)}
    for d in datas:
        for tau, rows in d['by_tau'].items():
            synthetic['by_tau'][tau].extend(rows)
    return summed_ratio(synthetic, label)


new_all = load('n16_all_p65536_d72', 'new-bitset')
new_one85 = load('n16_onepeak_p65536_d85', 'new-bitset')
new_one86 = load('n16_onepeak_p65536_d86', 'new-bitset')
new_tau2 = load('n16_tau2_p262144_d73', 'new-bitset')
old_524 = load('n12_p524288_d70', 'r013-old')
old_262_68 = load('n12_p262144_d68', 'r013-old')
old_262_69 = load('n12_p262144_d69', 'r013-old')

certificate = {
    'status': 'passed',
    'precision': 'Decimal(80)',
    'seed_check': {
        'new_bitset_source_formula': '(0xA022000000000000 + domain + n*1000000 + tau*10000 + replicate) mod 2^64',
        'r013_old_source_formula': '(0xA022000000000000 + seedDomain + n*1000000 + tau*10000 + replicate) mod 2^64',
        'finding': 'Verified independently from each file schema/source family. For these files the two displayed affine formulas are algebraically identical; domain field names and implementations differ.',
    },
    'inputs': [{k: v for k, v in d.items() if k in ('name','family','path','sha256','domain','dead_count','coverage')} | {
        'n': d['meta']['n'], 'particles': d['meta']['particles'], 'replicates': d['meta']['replicates'],
        'tau_lo': d['meta']['tau_lo'], 'tau_hi': d['meta']['tau_hi']}
        for d in (new_all,new_one85,new_one86,new_tau2,old_524,old_262_68,old_262_69)],
    'n16_doublepeak_P65536': summed_ratio(new_all, 'sum_tau mean_R; NOT uniform tau mixture'),
    'n16_tau2_P262144_R64': {'per_tau': per_tau_diagnostics(new_tau2)},
    'n16_onepeak_P65536_R256': {
        'domain85': {'per_tau': per_tau_diagnostics(new_one85)},
        'domain86': {'per_tau': per_tau_diagnostics(new_one86)},
        'same_P_pooled': {'per_tau': per_tau_diagnostics(combine_data := {'by_tau': defaultdict(list, {257: new_one85['by_tau'][257] + new_one86['by_tau'][257]})})},
    },
    'n12_P524288_R16': summed_ratio(old_524, 'P524288 only; no cross-P pooling'),
    'n12_P262144_R16_domains68_69': {
        'domain68': summed_ratio(old_262_68, 'P262144 domain68 only'),
        'domain69': summed_ratio(old_262_69, 'P262144 domain69 only'),
        'same_P_pooled': combine_same_p((old_262_68, old_262_69), 'P262144 domains68+69; same P only'),
    },
    'scope': 'All dead rows enter the raw Z/U means and sample variances as zeros. No n12 P524288/P262144 cross-particle pooling is performed. All Q and DeltaSE values are diagnostics only.',
}

ROOT.mkdir(parents=True, exist_ok=True)
(ROOT / 'diagnostics-certificate.json').write_text(json.dumps(certificate, indent=2) + '\n')
print(json.dumps({'status': certificate['status'], 'n16_Z2': certificate['n16_doublepeak_P65536']['Z_sum_tau_means'],
                  'n16_Q': certificate['n16_doublepeak_P65536']['Q_ratio'],
                  'n12_P524288_Q': certificate['n12_P524288_R16']['Q_ratio']}, indent=2))
