import pathlib,subprocess,sys,os,json,argparse,time,hashlib
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);z=p.parse_args();assert str(z.output.resolve()).startswith('/Volumes/data/CodexTemp/');z.output.mkdir(exist_ok=False)
b=pathlib.Path(__file__).resolve().parent;records=[]
def run(name,cmd,env=None):
 t=time.monotonic()
 with (z.output/(name+'.stdout')).open('w') as so,(z.output/(name+'.stderr')).open('w') as se:
  p=subprocess.Popen(cmd,stdout=so,stderr=se,env=env);rc=p.wait()
 try:os.kill(p.pid,0);absent=False
 except ProcessLookupError:absent=True
 records.append(dict(name=name,pid=p.pid,exit_code=rc,waited_reaped=True,absent_after=absent,wall_seconds=time.monotonic()-t));(z.output/'receipt.json').write_text(json.dumps(dict(python=sys.version,records=records),indent=2)+'\n');assert rc==0 and absent,name
run('ledger',[sys.executable,'-B',str(b/'eleven-prefix/audit_ledger.py')]);run('row',[sys.executable,'-B',str(b/'prefix-occupation/audit_row.py')])
inputs=[b.parent/'r036/finite-candidates'/('n%d-bfs-center-s%d.txt'%(n,s)) for n,s in [(8,2),(12,1)]]+sorted((b/'finite-words').glob('*.word'))
for i,w in enumerate(inputs):run('word%d'%i,[sys.executable,'-B',str(b/'prefix-occupation/analyze.py'),str(w),'--output',str(z.output/('word%d.json'%i))])
e=os.environ.copy();e['CLANG_MODULE_CACHE_PATH']='/Volumes/data/CodexTemp/cache/clang';e['TMPDIR']=str(z.output)
s=b.parent/'r037/large-side/large_side_orbit_dp.cpp';binary=z.output/'n11';run('compile',['clang++','-O3','-std=c++20','-Wall','-Wextra','-Werror','-pedantic',str(s),'-o',str(binary)],e)
run('n11-samples',[str(binary),'--mode','gate-n11-transpose','--out-dir',str(z.output/'gates'),'--source-sha',hashlib.sha256(s.read_bytes()).hexdigest()])
print('PASS',len(records),'checks; production enumeration not repeated')
