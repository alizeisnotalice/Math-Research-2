"""Exact rational feasibility-distance certificate for the reported decimal summaries."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent
INPUT=HERE.parent/'r025-genealogy'/'n24-genealogy-diagnostic.json'
x=json.loads(INPUT.read_text()); out={'scope':'Deterministic error floor for the recorded finite-decimal plug-in pair, not a confidence interval or a sampler-bug diagnosis','input_sha256':hashlib.sha256(INPUT.read_bytes()).hexdigest(),'runs':{}}
for name,run in x['runs'].items():
 t=run['terminal_weight_diagnostic']; b_hat=1/F(t['E_tau2_L'])
 p=list(map(F,t['complete_K2_L_weighted_distance_histogram']))
 assert p[0]==p[1]==0 and all(v>=0 for v in p)
 total=sum(p);p=[v/total for v in p]
 D=len(p)-1;g=lambda d:F(d,2*(d-1))
 defect=b_hat-sum(g(d)*p[d] for d in range(2,D+1))
 span=1-g(D); floor=defect/(1+span)
 threshold=F(15,100) if name=='local_strat' else F(18,100)
 assert floor>threshold
 out['runs'][name]={'histogram_renormalization':str(total),'defect':str(defect),'coefficient_span':str(span),'max_abs_normalizer_error_or_TV_error_lower_bound':str(floor),'strict_lower_bound_verified':str(threshold),'display_lower_bound':float(floor)}
(HERE/'plugin-error-floor.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
