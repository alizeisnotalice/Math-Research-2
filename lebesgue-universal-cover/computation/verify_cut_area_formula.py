"""Replay all dependent computations before accepting the fixed-angle area formula."""
from pathlib import Path
from fractions import Fraction as F
import runpy,io,contextlib,json,sys,hashlib
base=Path(__file__).parent

def replay(script,args):
 saved=sys.argv;buf=io.StringIO()
 try:
  sys.argv=[str(base/script),*args]
  with contextlib.redirect_stdout(buf):runpy.run_path(str(base/script),run_name='__main__')
 finally:sys.argv=saved
 return json.loads(buf.getvalue())
def mathematical(r):return {k:v for k,v in r.items() if k!='python'}
ah=replay('certify_ah_area.py',['2048'])
frozen_ah=json.loads((base/'022-ah-area-refined.json').read_text())
assert mathematical(ah)==mathematical(frozen_ah),'AH mathematical fields differ from actual replay'
cut=replay('certify_cut_area_formula.py',[])
frozen_cut=json.loads((base/'022-cut-area-formula.json').read_text())
assert mathematical(cut)==mathematical(frozen_cut),'Cut formula differs from actual replay'
lo,hi=map(F,cut['full_formula_interval']);assert lo<=hi<F('0.8440935944')
paths=['rational_trig.py','certify_ah_path.py','certify_ah_area.py','certify_cut_area_formula.py','022-ah-area-refined.json','022-cut-area-formula.json']
print(json.dumps({'python':sys.version,'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'AH_2048_reexecuted':True,'cut_formula_reexecuted':True,'all_mathematical_fields_match_frozen_outputs':True,'dependency_sha256':{p:hashlib.sha256((base/p).read_bytes()).hexdigest() for p in paths},'formula_interval':cut['full_formula_interval'],'strict_target':'0.8440935944','passed':True,'scope':'full arithmetic replay of the conditional cut-area formula, not geometric or universal-coverage verification'},indent=2))
