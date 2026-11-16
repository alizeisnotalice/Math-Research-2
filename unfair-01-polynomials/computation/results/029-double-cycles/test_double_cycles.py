"""Exercise the actual integer verifier with corrupted modular witnesses."""
import copy,json
from pathlib import Path
from verify_double_cycles import verify
def main():
 root=Path(__file__).resolve().parent;z=json.loads((root/'double-cycle-results.json').read_text());verify(z);passed=[]
 mutations=[('wrong_factor',lambda v:v['p3']['H_xp'][0]['ascending'].__setitem__(0,v['p3']['H_xp'][0]['ascending'][0]+1)),
            ('wrong_M',lambda v:v['p5'].__setitem__('M',5)),
            ('wrong_multiplicity',lambda v:v['p3']['H_xq'][0].__setitem__('multiplicity',2)),
            ('composite_ell',lambda v:v['p3'].__setitem__('ell',169))]
 for name,mutate in mutations:
  v=copy.deepcopy(z);mutate(v)
  try:verify(v)
  except (ValueError,AssertionError):passed.append(name)
  else:raise AssertionError(name)
 (root/'semantic-tests.json').write_text(json.dumps({'status':'production_integer_verifier_tampers_rejected','tests':passed},indent=2)+'\n');print(passed)
if __name__=='__main__':main()
