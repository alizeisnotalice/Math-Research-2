#!/usr/bin/env python3
"""R45: one R41-base q1 group[26,27] strictly-local slope certificate."""
import hashlib,json,os,subprocess,sys,time
from fractions import Fraction as F
from pathlib import Path
if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
p=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11');out=p/'round045/bracket-refine-q0';out.mkdir(parents=True,exist_ok=False)
inp=p/'round041/cross-group-probe/candidate_0036.json';base=p/'round008/joint-groups/batch2/candidate_0611.json';data=json.loads(base.read_text());target=json.loads(inp.read_text());assert len(data['pts'])==51 and len(target['pts'])==204
px=F(data['Px']);v=[0]*204;v[26]=1;v[27]=1;(out/'velocities.json').write_text(json.dumps(v)+'\n');s0=F("-139/172800000")*px; H=px/F(100000000)
enum=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-position-slab-impl/planar-centered-square-weak11/computation/position_slab_enumerator.py');sweep=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-operations-terra/r039_join_slope.py')
def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
cmd=[sys.executable,str(enum),str(inp),'--velocities',str(out/'velocities.json'),'--s0',str(s0),'--H',str(H),'--R','10','--output',str(out/'slab.json')];cmd2=[sys.executable,str(sweep),str(out/'slab.json'),'--output',str(out/'slope.json'),'--bits','80'];manifest={'scope':'R41 base q0 group original labels 26,27; physical s=s_normalized*base_px','normalized_s0':'-139/172800000','normalized_H':'1/100000000','physical_s0':str(s0),'physical_H':str(H),'input_sha256':sha(inp),'base_sha256':sha(base),'velocity_sha256':sha(out/'velocities.json'),'enum_sha256':sha(enum),'sweep_sha256':sha(sweep),'runner_sha256':sha(Path(__file__)),'commands':[cmd,cmd2]};(out/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for c in (cmd,cmd2):
 print('running',json.dumps(c),flush=True);r=subprocess.run(c,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True);print('exit',r.returncode,'stderr',r.stderr[-2000:],flush=True)
 if r.returncode:raise SystemExit(r.returncode)
print('complete',flush=True)
