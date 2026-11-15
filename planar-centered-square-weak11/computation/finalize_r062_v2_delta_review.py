#!/usr/bin/env python3
import hashlib,json,os,shutil
from pathlib import Path
assert os.path.ismount('/Volumes/data')
topic=Path(__file__).resolve().parents[1]
temp=Path('/Volumes/data/CodexTemp/T2/r062-independent-v2-audit')
dest=topic/'computation/results/round062-v2-independent-toys'
dest.mkdir(parents=True,exist_ok=True)
names=['toy-input.json','negative_incomplete.json','negative_complete.json','repaired_complete.json']
for name in names:shutil.copyfile(temp/name,dest/name)
result=topic/'computation/results/round062-v2-delta-independent.json'
data=json.loads(result.read_text());assert data['status']=='PASS_R62_V2_DELTA'
relative=['reviews/r062-v2-delta-independent-review.md','computation/review_r062_v2_delta.py',
          'computation/results/round062-v2-delta-independent.json','computation/finalize_r062_v2_delta_review.py']
relative+=['computation/results/round062-v2-independent-toys/'+n for n in names]
files=[]
for rel in relative:
    p=topic/rel;raw=p.read_bytes();assert len(raw)<10*1024*1024
    files.append(dict(source=str(p),destination=topic.name+'/'+rel,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw)))
out=dict(round=62,date='2026-11-15',status='V2 DELTA INDEPENDENT REVIEW PASS; V1 TERMINATION QUALIFICATION RESOLVED',
         scope='Tight derivative extrema, signed-curvature auxiliary LP, exact source/repair semantics; no optimized-candidate/full-D/new-C2 certificate.',
         files=files,total_payload_bytes=sum(f['bytes']for f in files),source_pins=data['source_pins'],
         actual_run=dict(python='3.13.18',sympy='1.14.0',exit_code=0,seconds=data['elapsed_seconds'],statistics=data['statistics'],
                         symbolic_identities=2,source_toy_runs=3),
         dependency='computation/review_r062_coupled_source.py from frozen v1 independent review supplies exact-Fraction reference functions.',
         git='No delegated commit/push; no MATLAB launched; root owns publication and final consolidated proof pin.')
p=topic/'computation/results/round062-v2-delta-independent-manifest.json';p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(path=str(p),files=len(files),bytes=out['total_payload_bytes'],sha256=hashlib.sha256(p.read_bytes()).hexdigest()),indent=2))
