#!/usr/bin/env python3
"""Mutation sensitivity checks for the exact R24 terminal-tree verifier.
No halo calculations; all changed certificate copies go to a data temp folder.
"""
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import time
from fractions import Fraction as F

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('passed_check', type=Path)
parser.add_argument('temp', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
start = time.process_time()
record = json.loads(args.passed_check.read_text())
command = record['actual_command']
checker_path = Path(command[1])
spec = importlib.util.spec_from_file_location('terminal_checker', checker_path)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)
paths = {command[i][2:]: Path(command[i + 1])
         for i in range(2, len(command), 2)}
paths.pop('output')
original = json.loads(paths['snapshot'].read_text())
args.temp.mkdir(parents=True, exist_ok=True)
mutations = {}

def reject(name, change):
    data = copy.deepcopy(original)
    change(data)
    p = args.temp / (name + '.json')
    p.write_text(json.dumps(data, separators=(',', ':')))
    altered = dict(paths, snapshot=p)
    try:
        checker.check(**altered)
    except AssertionError:
        mutations[name] = 'rejected'
        return
    raise AssertionError(name + ' unexpectedly accepted')

reject('missing_child_reference', lambda d: d['nodes'][0].update(children=[1, 1]))
reject('wrong_parent_backlink', lambda d: d['nodes'][1].update(parent=1))
reject('shifted_adjacent_boundary', lambda d: d['nodes'][1]['cell'].__setitem__(1, '-1/16'))
reject('incorrect_global_maximum', lambda d: d.update(upper=str(F(d['upper']) + 1)))
reject('incorrect_gap', lambda d: d.update(gap='0'))
reject('false_target_status', lambda d: d.update(status='target_gap'))
reject('outside_reference_point', lambda d: d.update(reference_point=['2', '0']))
reject('wrong_family_identity', lambda d: d.update(input_sha256='0' * 64))
result = {'status': 'passed', 'scope': 'malformed terminal-certificate rejection only; no halo geometry replay',
          'python': sys.version, 'platform': platform.platform(),
          'cpu_seconds': time.process_time() - start, 'mutations': mutations,
          'baseline_check_sha256': hashlib.sha256(args.passed_check.read_bytes()).hexdigest(),
          'checker_sha256': hashlib.sha256(checker_path.read_bytes()).hexdigest(),
          'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'actual_command': [sys.executable, str(Path(__file__).resolve()),
                             str(args.passed_check.resolve()),
                             str(args.temp.resolve()), str(args.output.resolve())]}
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'actual_command'}, indent=2))
