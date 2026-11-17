#!/usr/bin/env python3
"""Check expanded exact-Walsh receipts and report measured utilization."""
import argparse
import hashlib
import json
import re
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
base = Path(__file__).parent
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
records = []
for name, source, runner, per_l, seed0, workers in [
    ('remote-full-cpu-execution.json', 'cddl_fact22_lift_verify.cpp',
     'run_remote_full_cpu.py', 1024, 202611172000, 15),
    ('remote-qos-cpu-execution.json', 'cddl_lift_qos_verify.cpp',
     'run_remote_qos_cpu.py', 512, 202611173000, 15),
    ('remote-qos-30-execution.json', 'cddl_lift_qos_verify.cpp',
     'run_remote_qos_30.py', 256, 202611174000, 30),
]:
    receipt = json.loads((base/name).read_text())
    assert receipt['ssh_exit_code'] == 0
    assert receipt['runner_sha256'] == sha(base/runner)
    r = receipt['result']
    assert r['status'] == 'passed' and r['logical_cpus'] == 15
    assert r['source_sha256'] == sha(base/source)
    assert len(r['jobs']) == workers
    for i, j in enumerate(r['jobs']):
        s = j['result']
        assert j['exit_code'] == 0 and s['status'] == 'passed'
        assert s['seed'] == seed0+i
        assert s['sample_count'] == 5*per_l
        assert s['direct_coefficient_count'] == per_l*sum(1 << (14+2*l) for l in range(5))
        assert s['all_sample_predicates_verified']
        assert s['seed_walsh_all_abs_128']
        assert len(s['raw_result_sha256']) == 64
    lines = [line for s in r['whole_machine_top']['samples'] for line in s['lines']]
    idle = [float(re.search(r'([\d.]+)% idle', line)[1]) for line in lines]
    records.append({
        'receipt': name, 'receipt_sha256': sha(base/name), 'source_sha256': sha(base/source),
        'workers': workers, 'sample_count': workers*5*per_l,
        'direct_coefficient_count': sum(j['result']['direct_coefficient_count'] for j in r['jobs']),
        'wall_seconds': receipt['elapsed_seconds'],
        'task_cpu_peak_percent': max(s['percent_sum'] for s in r['cpu_samples']),
        'whole_machine_min_idle_percent': min(idle),
        'whole_machine_zero_idle_samples': sum(x == 0 for x in idle),
        'top_sample_count': len(idle),
    })
out = {'status': 'passed', 'batches': records,
       'scope': 'Finite integer Walsh audits; samples are not a new nonweak-normality certificate.',
       'source_sha256': sha(Path(__file__))}
a.out.write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out))
