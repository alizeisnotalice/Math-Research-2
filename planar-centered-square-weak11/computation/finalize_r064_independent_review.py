#!/usr/bin/env python3
import hashlib,json,os
from pathlib import Path
assert os.path.ismount('/Volumes/data')
topic=Path(__file__).resolve().parents[1]
a=json.loads((topic/'computation/results/round064-mass-variation-independent.json').read_text())
b=json.loads((topic/'computation/results/round064-continuation-competitor-independent.json').read_text())
assert a['status']=='PASS_R64_PROOF_GEOMETRY_AND_BEST_FINITE_STEP'
assert b['status']=='PASS_R64_CONTINUATION_ENDPOINT_AND_COMPETITOR_LEMMA'
relative=['reviews/r064-mass-polytope-independent-review.md','reviews/r064-competitor-envelope-independent-review.md',
          'computation/review_r064_mass_variation.py','computation/review_r064_continuation_competitor.py',
          'computation/results/round064-mass-variation-independent.json','computation/results/round064-continuation-competitor-independent.json',
          'computation/finalize_r064_independent_review.py']
files=[]
for rel in relative:
    p=topic/rel;raw=p.read_bytes();assert len(raw)<10*1024*1024
    files.append(dict(source=str(p),destination=topic.name+'/'+rel,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
proofs=[r for r in b['source_pins']if Path(r['path']).name in ('r064-source-state-mass-polytope-variation.md','r064-overlapping-competitor-envelope.md')]
for r in proofs:assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256']
out=dict(round=64,date='2026-11-15',status='INDEPENDENT PROOF/DOMAIN/SEED-BEST-ENDPOINT GEOMETRY REVIEW PASS',
         files=files,total_payload_bytes=sum(f['bytes']for f in files),final_reviewed_proofs=proofs,
         scope='General finite-source state and mass-polytope theorem plus competitor-envelope sufficient lemma reviewed. All72 saved masses independently feasible; all8 accepted continuation state/hash chains checked; complete geometry independently replayed ONLY originalseed, firstbest and finalendpoint.42 positive stored interval comparisons checked arithmetically, NOT all72 geometric certificates replayed. No gradient/local/global/C2 claim.',
         actual_runs=[dict(python='3.13.18',exit_code=0,seconds=a['elapsed_seconds'],geometry_cases=['seed','pilotbest']),
                      dict(python='3.13.18',exit_code=0,seconds=b['elapsed_seconds'],geometry_cases=['continuationendpoint'],
                           continuation_constraint_points=64,accepted_chain_steps=8,competitor_scalar_cases=4352)],
         final_gain_over_original_seed=b['independently_certified_final_gain_over_original_seed'],
         final_competitor_mass=b['final_L_star'],
         dependency='computation/review_r064_continuation_competitor.py imports the independent functions in computation/review_r064_mass_variation.py; neither imports producer.',
         git='No delegated commit/push, no MATLAB or optimizer launch; root owns archive. R62/R63 frozen files unchanged.')
p=topic/'computation/results/round064-independent-review-manifest.json';p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(path=str(p),files=len(files),bytes=out['total_payload_bytes'],sha256=hashlib.sha256(p.read_bytes()).hexdigest()),indent=2))
