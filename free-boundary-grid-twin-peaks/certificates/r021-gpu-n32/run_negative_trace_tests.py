#!/usr/bin/env python3
"""Mutate disposable copies; each independent verifier must reject it."""
import copy
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

base = Path(__file__).resolve().parent
trace = base / 'n32-one-final-p1-l2.jsonl.trace.json'
data = json.loads(trace.read_text())

cases = []
# Keep multi-megabyte altered traces out of the minimal archive.  The caller
# may set TASK_TMP; the default follows the project's data-volume convention.
tmp_root = Path(os.environ.get('TASK_TMP', '/Volumes/data/CodexTemp/T2'))
tmp_root.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix='a022-r021-negative-', dir=tmp_root) as td:
    work = Path(td)
    for name, mutate, verifier in (
        (
            'high_word',
            lambda x: x['rows'][-1]['state'][0][0].__setitem__(8, x['rows'][-1]['state'][0][0][8] ^ 1),
            base / 'verify_trace.py',
        ),
        (
            'truncated_prefix',
            lambda x: x['rows'][500]['prefix'][0].pop(),
            base / 'verify_trace.py',
        ),
        (
            'inverse_group',
            lambda x: x['rows'][-1]['inverse_m3_groups'].__setitem__(0, str(int(x['rows'][-1]['inverse_m3_groups'][0]) + 1)),
            base / 'verify_inverse.py',
        ),
    ):
        altered = copy.deepcopy(data)
        mutate(altered)
        out = work / f'negative-{name}.trace.json'
        out.write_text(json.dumps(altered))
        # The verifier derives the raw replicate file by removing .trace.json.
        raw = Path(str(out).replace('.trace.json', ''))
        raw.write_bytes((base / 'n32-one-final-p1-l2.jsonl').read_bytes())
        proc = subprocess.run([sys.executable, str(verifier), str(out)], text=True,
                              capture_output=True)
        assert proc.returncode != 0, (name, proc.stdout, proc.stderr)
        cases.append({'case': name, 'rejected': True, 'returncode': proc.returncode})

(base / 'negative-tests.json').write_text(json.dumps({'status': 'passed',
    'source_trace': trace.name, 'cases': cases}, indent=2) + '\n')
print(json.dumps({'status': 'passed', 'cases': cases}))
