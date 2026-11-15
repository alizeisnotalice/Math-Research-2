#!/usr/bin/env python3
"""Package bounded R62 independent audit and minimal complete toy receipts."""
import hashlib
import json
import os
from pathlib import Path
import shutil

assert os.path.ismount('/Volumes/data')
topic=Path(__file__).resolve().parents[1]
temp=Path('/Volumes/data/CodexTemp/T2/r062-independent-source-audit')
dest=topic/'computation/results/round062-independent-toys'
dest.mkdir(parents=True,exist_ok=True)
names=['invalid-toy-input.json','budget_incomplete.json','negative_early.json','negative_complete.json','repaired-toy.json']
for name in names:
    shutil.copyfile(temp/name,dest/name)
result=topic/'computation/results/round062-coupled-source-independent.json'
data=json.loads(result.read_text())
assert data['status']=='PASS_BOUNDED_INDEPENDENT_R62_AUDIT'
relative=['reviews/r062-coupled-source-independent-review.md','computation/review_r062_coupled_source.py',
          'computation/results/round062-coupled-source-independent.json','computation/finalize_r062_independent_review.py']
relative += ['computation/results/round062-independent-toys/'+n for n in names]
files=[]
for rel in relative:
    p=topic/rel;raw=p.read_bytes()
    assert len(raw)<10*1024*1024
    files.append(dict(source=str(p),destination=topic.name+'/'+rel,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
out=dict(round=62,date='2026-11-15',
         status='BOUNDED SOURCE CHECKER PASS; TERMINATION PROOF REQUIRES JOINT BITS/ATAN-TERMS QUALIFICATION',
         scope='Independent proof/code review and exact rational toy checks; no MATLAB launch, no optimized source/full-D/new-C2 certificate.',
         files=files,total_payload_bytes=sum(f['bytes']for f in files),
         author_sources_at_successful_run=data['source_pins'],actual_run=dict(python='3.13.18',sympy='1.14.0',exit_code=0,
         seconds=data['elapsed_seconds'],symbolic_identities=5,scalar_jet_points=60,rectangle_probability_points=108,
         source_coverage_repair_runs=4),git='No delegated commit or push; root owns integration.')
p=topic/'computation/results/round062-coupled-source-independent-manifest.json'
p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(path=str(p),files=len(files),bytes=out['total_payload_bytes'],sha256=hashlib.sha256(p.read_bytes()).hexdigest()),indent=2))
