from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
"""Replay balanced support certificates in both hands on every angle cell."""
from cross_common import *

raw=json.loads((D/'certificate-labels.json').read_text()); seeds=[tuple(x) for x in raw['seeds']];K=raw['K'];assert K==2048
expected={f'{h}-{b}-{s}' for h in ('r155','r175') for b in ('r151','r216') for s in (1,-1)}
assert set(raw)=={'K','seeds','rows'} and set(raw['rows'])==expected
assert len(seeds)==512 and len(set(seeds))==512
assert all(type(s) is int and s in (-1,1) and type(k) is int and abs(k)<=10**12 for s,k in seeds)
assert all(len(v)==K and all(len(z)==3 and len(set(z))==3 and all(type(j) is int and 0<=j<len(seeds) for j in z) for z in v) for v in raw['rows'].values())
pose155=json.loads((C/'free-p30-reoptimized-pose.json').read_text());pose175=json.loads((C/'seven-point-families-pose.json').read_text())['pose']
hosts={'r155':Host(pose155,False),'r175':Host(pose175,True)};cache={name:[h.row(s) for s in seeds] for name,h in hosts.items()};assert A('.294')**2 > A('.58')**2 - A('.25')
# Constant width one plus |x| < .58 gives |h'(angle)| < .294 a.e.
bound=A('.294')*pi/(2*K);out={};started=time.monotonic()
for key,rows in raw['rows'].items():
 hn,bod,sg=key.split('-',2);sg=int(sg);assert len(rows)==K;f=hp if bod=='r151' else hf;best=None
 for i,indices in enumerate(rows):
  rs=[cache[hn][j] for j in indices];weights=balance([r[0] for r in rs])
  for half in (0,1):
   a=pi*q(F(4*i+2*half+1,2*K));value=arb(0)
   for weight,(_,u,h) in zip(weights,rs):
    # Pull back R_a diag(1,sg) and cancel translation by exact balance.
    v=rot(-a,u);v=(v[0],sg*v[1]);value+=q(weight)*(f(v)-h)
   lb=(value-bound).lower();assert lb>0,(key,i,half,lb);best=lb if best is None or lb<best else best
 thresholds={'r155-r151':'.00154117','r155-r216':'.00149681','r175-r151':'.00145235','r175-r216':'.00141054'}
 assert best>A(thresholds[f'{hn}-{bod}'])
 out[key]=str(best);print(key,'lower',best,flush=True)
files=['free-p30-reoptimized-pose.json','seven-point-families-pose.json','finite-point-leaf-input.json','multi-host-points.json','deformed_reuleaux_support.py','certify_b9_padding_four_minimizer.py','certify_multi_host_points.py','finite-arc-reduction/certify_five_arc_band.py'];res=dict(lower_bounds=out,cells_per_hand=2*K, parent_cells=K,unit_radius_upper='.58',angular_loss=str(bound.upper()),precision_bits=ctx.prec,flint=flint.__version__,python=platform.python_version(),runtime_s=time.monotonic()-started,inputs_sha256={n:sha256((C/n).read_bytes()).hexdigest() for n in files},labels_sha256=sha256((D/'certificate-labels.json').read_bytes()).hexdigest());(D/'strict-result.json').write_text(json.dumps(res,indent=2)+'\n')
