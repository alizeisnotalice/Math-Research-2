"""Exact whole-hull core improvement and naive grid cost; not global exclusion."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,sys,runpy,io,contextlib,math
base=Path(__file__).parent
# Re-run all source-point and local arithmetic checks, not only a file hash.
with contextlib.redirect_stdout(io.StringIO()):
    old=runpy.run_path(str(base/'check_local_lower_box.py'))
r=old['r'];A=old['A'];eta=F(2,10**9)
assert F(99,70)**2>2
s=F(1,10000);delta=(F(99,70)+F(3,5))*s;r0=F(1,2)
alpha=(r0-delta)/(r0+eta);lower=alpha**2*A
assert lower>F(419,500)
pi_hi=F(3141592653589793239,10**18)
widths=[F(42,100)]*12+[F(48,100)]*2+[pi_hi/5]+[2*pi_hi/n for n in (7,9,11,19)]+[2*pi_hi]
counts=[math.ceil(w/(2*s)) for w in widths]
assert len(widths)==20
out={'python':sys.version,'parent_sha256':hashlib.sha256((base/'013-eight-body-certificate.json').read_bytes()).hexdigest(),'source_generator_reexecuted':True,'center':r['parameter_decimals_interpreted_exactly'],'dimension':20,'halfwidth':str(s),'delta':str(delta),'fixed_disk_inradius':str(r0),'alpha':str(alpha),'lower':str(lower),'lower_display_only':float(lower),'threshold':'419/500','volume_ratio_to_round015':2**20,'naive_tensor_grid_cell_counts':counts,'naive_tensor_grid_cell_product':str(math.prod(counts)),'scope':'one local box excluded; tensor grid count is for a naive uniform layout, not a lower bound on adaptive proof complexity','passed':True}
print(json.dumps(out,indent=2))
