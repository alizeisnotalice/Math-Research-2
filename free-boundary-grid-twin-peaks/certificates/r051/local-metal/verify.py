"""Replay the two finite length settings and display paired diagnostics."""
from pathlib import Path
import json
import math
import importlib.util

HERE = Path(__file__).resolve().parent
CERT = HERE.parents[1]
spec = importlib.util.spec_from_file_location('a022_metal_verify', CERT/'r047/n20-metal/verify.py')
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)

source = CERT/'r026/r026-merge-rb-calibration/source/tau2_merge_rb_calibration_r026.swift'
run = legacy.check(HERE/'n12.jsonl', HERE/'finished.json', HERE/'resource.txt', source)
a = json.loads((HERE/'earlier-l16-summary.json').read_text())
assert (a['n'],a['P'],a['R'],a['L']) == (12,131072,32,16)
assert a['audit_passed'] and a['domain'] != run['domain']
assert math.isclose(sum(a['distance_mRB']), 1, rel_tol=1e-12)
summary = {'long_run':run,
           'comparison':{'previous_P':a['P'],'previous_R':a['R'],'previous_L':a['L'],
                         'previous_domain':a['domain'],'previous_Q_mRB':a['Q_mRB'],
                         'previous_far_n_over_2':sum(a['distance_mRB'][7:]),
                         'same_random_seeds':False, 'same_replication_count':False},
           'scope':'Finite Monte Carlo sensitivity only, with no certified mixing or error bound'}
print(json.dumps(summary,indent=2))
