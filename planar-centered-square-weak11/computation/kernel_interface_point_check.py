"""Fresh exact-ball point obstruction; never a full C4 verification.

Usage: python -B kernel_interface_point_check.py INPUT_ROOT [CORE_ROOT]
Requires python-flint==0.9.0; input directory holds copied hl_core.py,
scaled.json, candidate.json, unit-origin-spline.json. No legacy file writes.
"""
import hashlib
import importlib.util
import json
import platform
from fractions import Fraction as F
from pathlib import Path
import sys
import flint

root=Path(sys.argv[1]).resolve()
core_root=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else root
spec=importlib.util.spec_from_file_location('kernel_audit_core',core_root/'hl_core.py')
core=importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)
core.ctx.prec=192
base=core.Spec(root/'scaled.json')
current=core.Spec(root/'candidate.json')
unit=core.Spec(root/'unit-origin-spline.json')
u,v=F(35,64),F(5,256)
def ak(s):
    return s.cA*s.APsi_point(u,v)+sum(s.AS_point(u,v))
def ball(x):
    return {'display':str(x),'lo':str(x.lower().fmpq()),'hi':str(x.upper().fmpq())}
ab,ac,ap=ak(base),ak(current),sum(unit.AS_point(u,v))
assert ap>0 and ac>0
assert (ac-(ab-core.A(F(1,200))*ap)).contains(0)
bad=ab-core.A(F(3,200))*ap
assert bad<0
cap=ac/ap
gain_cap=cap/core.A(2048)
assert gain_cap < core.A(F(1,200000))
mass=current.cA*current.a_*core.pw(current.RA,2-current.a_)/(2-current.a_)
mass+=core.A(current.h**2/2)*core.A(sum(current.coef(i,j) for i in range(-current.N,current.N+1) for j in range(-current.N,current.N+1)))
result={'scope':'one point necessary C4 constraint on the origin-spline ray only',
        'python':platform.python_version(),'python_flint':flint.__version__,'precision_bits':192,
        'inputs':{name:hashlib.sha256(((core_root if name=='hl_core.py' else root)/name).read_bytes()).hexdigest() for name in ['hl_core.py','scaled.json','candidate.json','unit-origin-spline.json']},
        'point':[str(u),str(v)],'AK_base':ball(ab),'AK_current':ball(ac),'AP_origin':ball(ap),
        'AK_base_minus_3over200_origin':ball(bad),'additional_epsilon_necessary_cap':ball(cap),
        'additional_halfmass_improvement_cap':ball(gain_cap),
        'additional_halfmass_improvement_strictly_below':'1/200000',
        'current_halfmass':ball(mass),
        'full_C1_C2_C4_rerun':False,'new_C2_upper_bound':False}
print(json.dumps(result,indent=2))
