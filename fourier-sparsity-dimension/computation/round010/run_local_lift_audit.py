#!/usr/bin/env python3
import argparse,hashlib,json,os,platform,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
a.work.mkdir(parents=True,exist_ok=False)
source=Path(__file__).with_name('cddl_fact22_lift_verify.cpp').resolve()
env=os.environ.copy();env.update(TMPDIR=str(a.work),TMP=str(a.work),TEMP=str(a.work))
compile_command=['clang++','-O3','-std=c++20','-fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang',str(source),'-o',str(a.work/'verify')]
c=subprocess.run(compile_command,capture_output=True,text=True,env=env);assert c.returncode==0,c.stderr
command=['nice','-n','10',str(a.work/'verify'),'--max-l','4','--samples','64','--seed','202611171000','--out',str(a.work/'result.json')]
beg=time.monotonic();r=subprocess.run(command,capture_output=True,text=True)
receipt={'host':platform.node(),'platform':platform.platform(),'compiler':subprocess.check_output(['clang++','--version'],text=True),
 'compile_command':compile_command,'command':command,'exit_code':r.returncode,'stdout':r.stdout,'stderr':r.stderr,
 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'binary_sha256':hashlib.sha256((a.work/'verify').read_bytes()).hexdigest(),
 'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-beg}
if r.returncode==0:receipt['result']=json.loads((a.work/'result.json').read_text())
a.out.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'exit_code':r.returncode,'elapsed_seconds':receipt['elapsed_seconds']}));raise SystemExit(r.returncode)
