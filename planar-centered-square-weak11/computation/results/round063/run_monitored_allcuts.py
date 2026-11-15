import subprocess,pathlib,os,time,json,sys,hashlib,concurrent.futures
b=pathlib.Path(__file__).resolve().parent
assert os.path.ismount('/Volumes/data')
o=b/'allcuts-batch';o.mkdir(exist_ok=True)
env=os.environ.copy()
for k in ['TMPDIR','TMP','TEMP','XDG_CACHE_HOME']:env[k]='/Volumes/data/CodexTemp/T2/frontier063'
pathlib.Path(env['TMPDIR']).mkdir(parents=True,exist_ok=True)
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
def run(s):
 cmd=[sys.executable,'-B',str(b/'exact_source_frontier.py'),'--input',str(b/'structural-inputs.json'),'--stratum',str(s),'--out',str(o/f'stratum-{s:02}.json')]
 r={'stratum':s,'command':cmd,'source_sha256':hashlib.sha256((b/'exact_source_frontier.py').read_bytes()).hexdigest(),'samples':[]};t=time.monotonic()
 with (o/f'stratum-{s:02}.log').open('w') as f:
  p=subprocess.Popen(cmd,env=env,cwd=b,stdout=f,stderr=subprocess.STDOUT);r['pid']=p.pid
  while p.poll() is None:
   q=subprocess.run(['ps','-p',str(p.pid),'-o','pcpu=,rss='],capture_output=True,text=True).stdout.strip();r['samples'].append({'elapsed':time.monotonic()-t,'ps':q})
   if (q and int(q.split()[-1])>768*1024) or time.monotonic()-t>1800:
    r['terminated_budget']=True;p.terminate()
    try:p.wait(timeout=5)
    except subprocess.TimeoutExpired:p.kill()
   time.sleep(.2)
  r['exit_status']=p.wait()
 r['wall_seconds']=time.monotonic()-t;r['owned_pid_absent']=not subprocess.run(['ps','-p',str(p.pid),'-o','pid='],capture_output=True,text=True).stdout.strip()
 out=o/f'stratum-{s:02}.json'
 if out.exists():r['output_sha256']=hashlib.sha256(out.read_bytes()).hexdigest()
 (o/f'receipt-{s:02}.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ['samples','command']}),flush=True);return r
t=time.monotonic()
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:rs=list(ex.map(run,range(12)))
(o/'batch-receipt.json').write_text(json.dumps({'workers':12,'wall_seconds':time.monotonic()-t,'jobs':rs},indent=2)+'\n')
assert all(r['exit_status']==0 and r['owned_pid_absent'] for r in rs)
