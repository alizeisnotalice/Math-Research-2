"""Verify submitted supports; this does not claim full layer coverage."""
import json
from pathlib import Path
from diagonal_domination import verify_certificate
from check_layer19_exceptions import check
base=Path(__file__).with_name('results')
cases=json.loads((base/'019-candidate-certificates.json').read_text())['cases']
key=lambda c:(tuple(c['I']),tuple(c['J']))
assert len(set(map(key,cases)))==len(cases)
special={key(c) for c in check()['cases']};seen=set();direct=0
for c in cases:
 if c['certificate'] is not None:
  verify_certificate(c['I'],c['J'],c['certificate']);direct+=1
 else:
  assert key(c) in special;seen.add(key(c))
assert seen==special
print(json.dumps({'all_submitted_supports_verified':True,'submitted':len(cases),'direct':direct,'special':len(seen),'full_layer_coverage_claimed':False},indent=2))
