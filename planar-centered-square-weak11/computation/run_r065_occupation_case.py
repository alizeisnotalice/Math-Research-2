#!/usr/bin/env python3
"""Reproduce ONE archived R65 case in a fresh data-volume directory.

This common wrapper is prepared after the five archived actual root runs;
those runs used the individual root wrapper hashes in actual-cases.json.
It does not assert those runs used this new wrapper.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,signal,subprocess,sys,time
BASE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert os.path.ismount('/Volumes/data')
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',required=True);ap.add_argument('--output-root',required=True)
    ap.add_argument('--matlab',default='/Applications/MATLAB_R2026b.app/bin/matlab')
    ap.add_argument('--max-seconds',type=float,default=240)
    ap.add_argument('--verify-source',action='store_true');ap.add_argument('--export-sparse',action='store_true');args=ap.parse_args()
    data=json.loads((BASE/'inputs/r065-tensor-occupation/actual-cases.json').read_text())
    spec=next(c for c in data['cases'] if c['id']==args.case)
    out=Path(args.output_root).resolve()
    assert str(out).startswith('/Volumes/data/') and not out.exists(),'Use a fresh data-volume directory'
    out.mkdir(parents=True)
    source=BASE/spec['matlab_source'];checker=BASE/'r065_tensor_source_verify.py'
    grid=BASE/'inputs/r065-tensor-occupation'/spec['grid_filename']
    assert sha(source)==spec['matlab_source_sha256'] and sha(grid)==spec['grid_sha256']
    if args.verify_source:assert 'source_check_limits' in spec and sha(checker)==spec['source_checker_sha256']
    if args.export_sparse:assert spec['function'].endswith('_v3')
    shutil.copyfile(source,out/source.name);shutil.copyfile(checker,out/checker.name);shutil.copyfile(grid,out/grid.name)
    (out/'input-case.json').write_text(json.dumps(spec,indent=2)+'\n')
    cache=Path('/Volumes/data/CodexTemp/T2/r065-kernel');cache.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    for key in ['TMPDIR','TMP','TEMP','XDG_CACHE_HOME']:
        if not env.get(key,'').startswith('/Volumes/data/'):env[key]=str(cache)
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[key]='1'
    # Use direct argv (no shell), with MATLAB string quoting for paths.
    quote=lambda x:"'"+str(x).replace("'","''")+"'"
    params=[format(spec['h'],'.17g'),quote(out/grid.name),str(spec['subdiv']),quote(out)]
    if spec['run_mode'] is not None:params.append(quote(spec['run_mode']))
    if spec['function'].endswith('_v3'):params.append('true'if args.export_sparse else'false')
    expr=spec['function']+'('+','.join(params)+')'
    cmd=[args.matlab,'-singleCompThread','-batch',expr]
    rec=dict(command=cmd,cwd=str(out),source_sha256=sha(source),grid_sha256=sha(grid),wrapper_sha256=sha(Path(__file__)),samples=[],thread_scope='Requested -singleCompThread and BLAS/OMP=1. Recorded CPU percent is not a proof every process phase uses one thread.')
    started=time.monotonic()
    with (out/'stdout.log').open('w') as log:
        proc=subprocess.Popen(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        rec['pid']=proc.pid
        while proc.poll() is None:
            ps=subprocess.run(['ps','-p',str(proc.pid),'-o','pcpu=,rss='],capture_output=True,text=True)
            rec['samples'].append(dict(elapsed=time.monotonic()-started,ps=ps.stdout.strip()))
            if time.monotonic()-started>args.max_seconds:
                rec['terminated_timeout']=True
                os.killpg(proc.pid,signal.SIGTERM)
                try:proc.wait(timeout=10)
                except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL)
                break
            time.sleep(1)
        rec['exit_status']=proc.wait()
    rec['wall_seconds']=time.monotonic()-started
    rec['owned_pid_absent']=not subprocess.run(['ps','-p',str(proc.pid),'-o','pid='],capture_output=True,text=True).stdout.strip()
    candidate=out/spec.get('candidate_filename','occupation-tensor-jump-lp.json')
    if candidate.exists():rec['output_sha256']=sha(candidate)
    (out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps({k:v for k,v in rec.items() if k!='samples'},indent=2))
    if rec['exit_status']!=0:return rec['exit_status'] or 1
    if args.verify_source:
        lim=spec['source_check_limits']
        vcmd=[sys.executable,str(out/checker.name),'--input',str(candidate),'--output',str(out/'source-certificate.json'),'--initial-subdiv',str(lim['initial_subdiv']),'--max-depth',str(lim['max_depth']),'--max-boxes',str(lim['max_boxes']),'--max-seconds',str(lim['max_seconds'])]
        result=subprocess.run(vcmd,cwd=out,env=env,text=True,capture_output=True)
        (out/'source-verify.log').write_text(result.stdout+result.stderr)
        print(result.stdout,end='')
        if result.returncode:return result.returncode
        certificate=json.loads((out/'source-certificate.json').read_text())
        return 0 if certificate['status']=='PASS_SOURCE_NONNEGATIVE' else 2
    return 0
if __name__=='__main__':raise SystemExit(main())
