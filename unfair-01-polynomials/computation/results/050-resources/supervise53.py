"""Wait/reap exact semantic runner and atomically preserve terminal status."""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path('/Users/zhengzhihao/CodexRemoteMath/a30-semantic53')
def atom(p,x):
    tmp=p.with_suffix('.partial')
    with tmp.open('w') as f:json.dump(x,f,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(p)
def main():
    command=[sys.executable,str(ROOT/'code/runner.py'),'--root',str(ROOT),'--total','53','--workers','13','--seconds','14400']
    with (ROOT/'compute53.log').open('a') as f:
        child=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT)
        atom(ROOT/'launch53.json',{'supervisor_pid':os.getpid(),'pid':child.pid,'command':command,'runner_sha256':hashlib.sha256((ROOT/'code/runner.py').read_bytes()).hexdigest(),'started_unix':time.time(),'scope':'Only total53 independently complete slice; original53/54 production interrupted.'})
        rc=child.wait()
    q=json.loads((ROOT/'progress53.json').read_text())
    atom(ROOT/'exit53.json',{'returncode':rc,'child_waited':True,'verified':q['verified'],'remaining':q['remaining'],'status':q['status'],'pid':child.pid,'ended_unix':time.time()})
if __name__=='__main__':main()
