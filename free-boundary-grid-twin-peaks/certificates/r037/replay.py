#!/usr/bin/env python3
"""Bounded exact checks, no production scan. Outputs must be on data volume."""
import argparse,hashlib,json,os,pathlib,subprocess,sys,time
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);z=p.parse_args()
assert str(z.output.resolve()).startswith('/Volumes/data/CodexTemp/')
z.output.mkdir(parents=True,exist_ok=False);base=pathlib.Path(__file__).resolve().parent
records=[]
def run(name,cmd,env=None):
 start=time.monotonic()
 with (z.output/(name+'.stdout')).open('w') as so,(z.output/(name+'.stderr')).open('w') as se:
  c=subprocess.Popen(cmd,stdout=so,stderr=se,env=env);rc=c.wait()
 try:os.kill(c.pid,0);absent=False
 except ProcessLookupError:absent=True
 records.append(dict(name=name,command=cmd,pid=c.pid,exit_code=rc,waited_reaped=True,absent_after=absent,wall_seconds=time.monotonic()-start))
 (z.output/'receipt.json').write_text(json.dumps(dict(python=sys.version,checks=records),indent=2)+'\n')
 if rc or not absent:raise SystemExit('check failed: '+name)
e=os.environ.copy();e['CLANG_MODULE_CACHE_PATH']='/Volumes/data/CodexTemp/cache/clang';e['TMPDIR']=str(z.output)
run('low-small',[sys.executable,'-B',str(base/'low-small/audit_low_small_orientation_n3.py')])
run('merge-depth',[sys.executable,'-B',str(base/'merge-depth/audit_path_merge_depth.py')])
for name,src in [('rollback',base/'undo-search/audit_rollback.cpp'),('enumerator',base/'large-side/audit_enumerator.cpp')]:
 binary=z.output/name
 run(name+'-compile',['clang++','-std=c++20','-O2','-Wall','-Wextra','-Werror',str(src),'-o',str(binary)],e)
 run(name,[str(binary)])
e['REPLAY_OUT']=str(z.output/'undo-gate');e['STEPS']='10000'
run('undo-small-gate',['bash',str(base/'undo-search/replay_small_gate.sh')],e)
print('PASS',len(records),'checks')
