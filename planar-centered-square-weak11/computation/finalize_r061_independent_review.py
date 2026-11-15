#!/usr/bin/env python3
"""Package this completed independent review without touching author files."""
import hashlib
import json
import os
from pathlib import Path
import shutil

assert os.path.ismount('/Volumes/data')
topic = Path(__file__).resolve().parents[1]
result = topic / 'computation/results/round061-source-union-independent.json'
data = json.loads(result.read_text())
assert data['status'] == 'PASS'
receipt_source = Path('/Volumes/data/CodexTemp/T2/r061-independent-audit/20261115T103253Z-39bd4809/run.json')
receipt = json.loads(receipt_source.read_text())
assert receipt['state'] == 'succeeded' and receipt['exit_code'] == receipt['child_returncode'] == 0
assert receipt['script_sha256'] == data['script_sha256']
receipt_dest = topic / 'computation/results/round061-source-union-independent-execution.json'
shutil.copyfile(receipt_source, receipt_dest)

relative = [
    'reviews/r061-source-union-independent-review.md',
    'computation/review_r061_source_union.py',
    'computation/results/round061-source-union-independent.json',
    'computation/results/round061-source-union-independent-execution.json',
    'computation/finalize_r061_independent_review.py',
]
files = []
for rel in relative:
    path = topic / rel
    raw = path.read_bytes()
    assert len(raw) < 10 * 1024 * 1024
    files.append(dict(source=str(path), destination=topic.name + '/' + rel,
                      bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
manifest = dict(
    round=61, date='2026-11-15', status='INDEPENDENT REVIEW PASS WITH STATED FINITE-INPUT SCOPE',
    scope='Proof reviewed; every saved geometry/cut/flow audited independently. No new general C2 bound.',
    files=files, total_payload_bytes=sum(r['bytes'] for r in files),
    author_manifest_path=data['manifest_path'], author_manifest_sha256=data['manifest_sha256'],
    independent_actual_run=dict(python=receipt['runtime']['version'], exit_code=receipt['exit_code'],
                                script_seconds=data['wall_seconds'], receipt_seconds=receipt['elapsed_seconds'],
                                statistics=data['statistics'], edge_inputs=len(data['edge_tests']),
                                edge_runs=sum(len(t['runs']) for t in data['edge_tests'])),
    inherited_source_status='R53 author theorem and independent Borel/size-inflation proof read and hash-checked; R53 certificate NOT rerun.',
    git='No delegated commit or push; root owns integration and publication.',
)
output = topic / 'computation/results/round061-source-union-independent-manifest.json'
output.write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps(dict(path=str(output), files=len(files), bytes=manifest['total_payload_bytes'],
                     sha256=hashlib.sha256(output.read_bytes()).hexdigest(), status=manifest['status']), indent=2))
