#!/usr/bin/env python3
"""Recheck the frozen r013 GPU batches and source/seed bookkeeping."""
from pathlib import Path
import hashlib
import json
import statistics
from audit_tau_smc_aggregates import aggregate, read, ROOT

BASE = ROOT/'certificates/r013-gpu'


def metal(path):
    return path.read_text().split('let metal = """', 1)[1].split('"""', 1)[0]


if __name__ == '__main__':
    for x in json.loads((BASE/'source-index.json').read_text()):
        assert hashlib.sha256((BASE/x['path']).read_bytes()).hexdigest() == x['sha256']
    original = metal(BASE/'audit/tau_smc_r013.swift')
    assert metal(BASE/'stratified/tau_smc_pool.swift') == original
    assert metal(BASE/'calibration/production_r013.swift') == original
    audit = json.loads((BASE/'audit/audit-results.json').read_text())
    assert audit['cases'] == 32 and audit['failures'] == 0
    rows = []
    seen = set()
    for path, domain, base in [
        ('calibration/production/n7-p65536-r16-strat-domain6300.jsonl', 0x6300000000000000, 0xA022000000000000),
        ('multinomial/raw/n8-p131072-r32.jsonl', 0x5100000000000000, 0),
        ('stratified/n8-p131072-r32-strat.jsonl', 0x7600000000000000, 0xA022000000000000),
    ]:
        f = BASE/path
        m, groups, raw = read(f)
        assert raw[-1]['kind'] == 'run_end'
        for t, reps in groups.items():
            for r, x in reps.items():
                seed = (base + domain + m['n']*1000000+t*10000+r) % (1 << 64)
                assert int(x['seed_hex'], 16) == seed and seed not in seen
                seen.add(seed)
        result, _ = aggregate(f)
        rows.append(result)
    telemetry = [json.loads(s) for s in (BASE/'stratified/n8-p131072-r32-strat-telemetry.jsonl').read_text().splitlines()]
    assert telemetry[-1]['kind'] == 'exit' and telemetry[-1]['returncode'] == 0
    v = [x['device_utilization_pct'] for x in telemetry if x['kind'] == 'sample' and x['device_utilization_pct'] is not None]
    print(json.dumps({'status':'passed', 'rows':rows, 'distinct_recorded_seeds':len(seen),
                      'whole_device_nonexclusive_telemetry':{'samples':len(v), 'mean_percent':statistics.mean(v),
                         'min_percent':min(v), 'max_percent':max(v), 'samples_at_100':v.count(100)}}, indent=2))
