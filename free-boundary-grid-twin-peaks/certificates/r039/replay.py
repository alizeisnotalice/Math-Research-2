import pathlib,argparse,subprocess,os,json,time,sys
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);z=p.parse_args();assert str(z.output.resolve()).startswith('/Volumes/data/CodexTemp/');z.output.mkdir(exist_ok=False)
b=pathlib.Path(__file__).resolve().parent;records=[]
def run(name,cmd,env=None):
 t=time.monotonic()
 with (z.output/(name+'.stdout')).open('w') as so,(z.output/(name+'.stderr')).open('w') as se:
  c=subprocess.Popen(cmd,stdout=so,stderr=se,env=env);rc=c.wait()
 try:os.kill(c.pid,0);absent=False
 except ProcessLookupError:absent=True
 records.append(dict(name=name,pid=c.pid,exit_code=rc,waited_reaped=True,absent_after=absent,wall_seconds=time.monotonic()-t));(z.output/'receipt.json').write_text(json.dumps(dict(python=sys.version,records=records),indent=2)+'\n');assert rc==0 and absent,name
run('mixed',[sys.executable,'-B',str(b/'mixed-mechanisms/audit_mixed_terms.py')])
run('parents',[sys.executable,'-B',str(b/'chain-obstruction/audit_eligible_parents.py')])
e=os.environ.copy();e['CLANG_MODULE_CACHE_PATH']='/Volumes/data/CodexTemp/cache/clang';e['TMPDIR']=str(z.output)
for tag,src in [('old',b.parent/'r037/undo-search/psi_search_undo.cpp'),('new',b/'warm-search/psi_search_warm.cpp')]:run('compile-'+tag,['clang++','-std=c++20','-O2','-Wall','-Wextra','-Werror',str(src),'-o',str(z.output/tag)],e)
run('warm-gates',[sys.executable,'-B',str(b/'warm-search/gate.py'),'--output',str(z.output/'gates'),'--old',str(z.output/'old'),'--new',str(z.output/'new')])
for i,w in enumerate(sorted((b/'finite-inputs').glob('*.word'))):run('word%d'%i,[sys.executable,'-B',str(b.parent/'r036/finite-candidates/verify_psi.py'),str(w)])
print('PASS',len(records),'root checks; warm gate includes18 processes')
