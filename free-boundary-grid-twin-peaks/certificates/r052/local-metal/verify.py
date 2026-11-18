"""Verify the n=16 Metal run and compare it with the earlier finite diagnostic."""
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import math

D=Path(__file__).resolve().parent
CERT=D.parents[1]
spec=importlib.util.spec_from_file_location('a022_calibration',CERT/'r047/n20-metal/verify.py')
legacy=importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)
source=CERT/'r026/r026-merge-rb-calibration/source/tau2_merge_rb_calibration_r026.swift'
new=legacy.check(D/'n16.jsonl',D/'finished.json',D/'resource.txt',source)
control=D/'control'
short=legacy.check(control/'n16.jsonl',control/'finished.json',control/'resource.txt',source)
assert (new['n'],new['P'],new['R'],new['domain']) == (short['n'],short['P'],short['R'],short['domain'])
assert new['seed_hex']==short['seed_hex'] and (short['L'],new['L'])==(16,128)
earlier=CERT/'r028/r028-merge-rb-n16'
ledger=earlier/'reports/n16-production-audit.json'
old=json.loads(ledger.read_text())
assert old['audit_passed'] and (old['n'],old['P'],old['R'],old['L'])==(16,262144,64,16)
assert math.isclose(sum(old['distance_mRB']),1,rel_tol=1e-12)
raw=earlier/'raw/n16-p262144-r64-l16.jsonl'
assert sha256(raw.read_bytes()).hexdigest()==old['raw_sha256']
assert (earlier/'source/tau2_merge_rb_calibration_r026.swift').read_bytes()==source.read_bytes()
rows=[json.loads(x) for x in raw.read_text().splitlines()]
r,rb=rows[1:-1:2],rows[2:-1:2]
assert len(r)==64==len(rb)
base=max(x['logZ'] for x in r)
weights=[math.exp(x['logZ']-base) for x in r]
masses=[weight*x['mean_m'] for weight,x in zip(weights,rb,strict=True)]
old_tail=sum(weight*sum(x['mean_m_hist'][9:]) for weight,x in zip(weights,rb,strict=True))/sum(masses)
assert math.isclose(old_tail,sum(old['distance_mRB'][9:]),rel_tol=1e-11)
assert new['domain']!=old['domain']
print(json.dumps({'long_run':new,'matched_short_run':short,
      'long_minus_matched_short':{'far_n_over_2':new['far_n_over_2']-short['far_n_over_2'],
                                  'Q_mRB':new['Q_mRB']-short['Q_mRB']},
      'earlier':{'P':old['P'],'R':old['R'],'L':old['L'],
      'seed_domain':old['domain'],'Q_mRB':old['Q_mRB'],'far_n_over_2':old_tail,
      'replicate_weight_ess':sum(masses)**2/sum(x*x for x in masses),
      'normalized_rep_weight_max':max(masses)/sum(masses)},
      'scope':'Matched settings for new runs except the number of moves; finite sensitivity only, no mixing bound'},
      indent=2))
