from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,os,json,time,hashlib

ROOT=Path('/Volumes/data/CodexProjects/T2');BASE=Path(__file__).parent
assert os.path.ismount('/Volumes/data')
SCRIPT=BASE/'theta_competitor_budget.py'
assert hashlib.sha256(SCRIPT.read_bytes()).hexdigest()=='211d33e21e47e6e14290e4ddbc3b3a36e247131de1bc53114091b7766d0fe8b6'
CACHE=Path('/Volumes/data/CodexTemp/T2/frontier066');CACHE.mkdir(parents=True,exist_ok=True)
ENV=os.environ.copy()
for k in ('TMPDIR','TMP','TEMP','XDG_CACHE_HOME'):
    if not ENV.get(k,'').startswith('/Volumes/data/'):ENV[k]=str(CACHE)
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):ENV[k]='1'

def run(i):
    cmd=['/Users/zzh/.venvs/main/bin/python3','-B',str(SCRIPT),'--module',str(ROOT/'.local/compute/r064-lp-pilot/mass_polytope_variation.py'),'--r63-input',str(ROOT/'.local/compute/frontier063/structural-inputs.json'),'--stratum',str(i),'--out',str(BASE/f'stratum-{i:02}.json')]
    t=time.monotonic();samples=[]
    with (BASE/f'stratum-{i:02}.stdout').open('w') as log:
        p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=ENV)
        while p.poll() is None:
            samples.append(subprocess.run(['ps','-p',str(p.pid),'-o','pcpu=,rss='],capture_output=True,text=True).stdout.strip())
            if time.monotonic()-t>180:
                p.terminate()
                try:p.wait(timeout=5)
                except subprocess.TimeoutExpired:p.kill()
                break
            time.sleep(.05)
        ec=p.wait()
    receipt={'stratum':i,'command':cmd,'pid':p.pid,'exit_status':ec,'wall_seconds':time.monotonic()-t,'owned_pid_absent':not subprocess.run(['ps','-p',str(p.pid),'-o','pid='],capture_output=True,text=True).stdout.strip(),'samples':samples}
    (BASE/f'stratum-{i:02}.receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    assert ec==0,(i,(BASE/f'stratum-{i:02}.stdout').read_text())
    data=json.loads((BASE/f'stratum-{i:02}.json').read_text());assert data['status']=='PASS'
    return {'stratum':i,'blocks':data['evaluated_proper_blocks'],'cases':len(data['cases']),'receipt':receipt}

t=time.monotonic()
with ThreadPoolExecutor(max_workers=11) as pool:rows=list(pool.map(run,range(12)))
result={'status':'PASS','root_parallel_jobs':11,'reason':'one concurrent root MATLAB batch; twelve fixed strata queued','wall_seconds':time.monotonic()-t,'total_cases':sum(r['cases'] for r in rows),'total_blocks':sum(r['blocks'] for r in rows),'rows':rows}
(BASE/'batch-receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
