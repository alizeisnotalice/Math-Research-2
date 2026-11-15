#!/usr/bin/env python3
"""Fixed t=1/2, physical z=Px/10^6; R39 one-axis audit, not a 2D result."""
import hashlib,json,os,subprocess,sys
from fractions import Fraction as F
from pathlib import Path
if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
here=Path(__file__).resolve().parent;slice=here/'slice';out=here/'one_axis';out.mkdir(exist_ok=False);m=json.loads((slice/'manifest.json').read_text());px=F(m['px']);enum=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-position-slab-impl/planar-centered-square-weak11/computation/position_slab_enumerator.py');sweep=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-operations-terra/r039_join_slope.py');sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();cmd=[sys.executable,str(enum),str(slice/'midmass_t_half.json'),'--velocities',str(slice/'velocity_x.json'),'--s0',str(px/F(1000000)),'--H',str(px/F(100000000)),'--R','10','--output',str(out/'slab.json')];cmd2=[sys.executable,str(sweep),str(out/'slab.json'),'--output',str(out/'slope.json'),'--bits','80'];(out/'run_manifest.json').write_text(json.dumps({'scope':'fixed t=1/2 only; no common 2D stratum claim','normalized_z0':'1/1000000','normalized_H':'1/100000000','slice_manifest_sha256':sha(slice/'manifest.json'),'enumerator_sha256':sha(enum),'sweep_sha256':sha(sweep),'runner_sha256':sha(__file__),'commands':[cmd,cmd2]},indent=2)+'\n')
for c in (cmd,cmd2):
 print('running',json.dumps(c),flush=True);r=subprocess.run(c,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True);print('exit',r.returncode,'stderr',r.stderr[-2000:],flush=True)
 if r.returncode:raise SystemExit(r.returncode)
print('complete',flush=True)
