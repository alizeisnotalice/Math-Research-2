#!/usr/bin/env python3
import hashlib,json,os,shutil
from pathlib import Path
assert os.path.ismount('/Volumes/data')
topic=Path(__file__).resolve().parents[1]
temp=Path('/Volumes/data/CodexTemp/T2/r062-independent-v3-audit')
dest=topic/'computation/results/round062-v3-independent-toys';dest.mkdir(parents=True,exist_ok=True)
names=['toy-input.json','negative_incomplete.json','negative_complete.json','repaired_complete.json']
for name in names:shutil.copyfile(temp/name,dest/name)
data=json.loads((topic/'computation/results/round062-v3-jump-independent.json').read_text())
assert data['status']=='PASS_R62_V3_JUMP_DELTA'
relative=['reviews/r062-v3-jump-independent-review.md','computation/review_r062_v3_jump_delta.py',
          'computation/results/round062-v3-jump-independent.json','computation/finalize_r062_v3_jump_review.py']
relative+=['computation/results/round062-v3-independent-toys/'+name for name in names]
files=[]
for rel in relative:
    p=topic/rel;raw=p.read_bytes();assert len(raw)<10*1024*1024
    files.append(dict(source=str(p),destination=topic.name+'/'+rel,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw)))
out=dict(round=62,date='2026-11-15',status='V3 JUMP DELTA INDEPENDENT REVIEW PASS',
         scope='Telescoped signed physical density jumps, exterior zero terms, LP sign split and exact checker delta; no MATLAB launch, optimized-candidate/full-D/new-C2 certificate.',
         files=files,total_payload_bytes=sum(f['bytes']for f in files),source_pins=data['source_pins'],
         reviewed_proof=dict(path='/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11/proofs/r062-jump-telescoped-source-curvature.md',sha256='49eccfb093f8eb6912e7c64805baccfe811f32d007a3f17b139815b34693eee1'),
         actual_run=dict(python='3.13.18',sympy='1.14.0',exit_code=0,seconds=data['elapsed_seconds'],symbolic_identities=2,
                         fullcell_curvature_points=108,boundary_checks=22,source_coverage_repair_runs=3),
         dependency='Frozen computation/review_r062_coupled_source.py supplies independent exact-Fraction arctan reference.',
         git='No delegated commit/push; root owns final archive and consolidated source pin.')
p=topic/'computation/results/round062-v3-jump-independent-manifest.json';p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(path=str(p),files=len(files),bytes=out['total_payload_bytes'],sha256=hashlib.sha256(p.read_bytes()).hexdigest()),indent=2))
