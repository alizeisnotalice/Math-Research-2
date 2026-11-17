#!/usr/bin/env python3
"""Strict, fresh 47/48 bridge.  Preparation does not invoke this script."""
from __future__ import annotations
import argparse,base64,gzip,hashlib,json,os,subprocess,time
from pathlib import Path
from support48_prepare import FROZEN_SHA,sha,validate

HERE=Path(__file__).resolve().parent
SSH_OPTIONS=["-o","HostName=fe80::8c00:aeff:feee:2fe%%anpi1","-o","HostKeyAlias=169.254.63.75","-o","ControlMaster=auto","-o","ControlPersist=600","-o","ControlPath=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/macbook-anpi1-20261117.sock","-o","BatchMode=yes","-o","StrictHostKeyChecking=yes","-o","ConnectTimeout=10"]
REMOTE_COMMAND="nice -n 10 env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/zhengzhihao/.venvs/math/bin/python -u -"

def sync(f): f.flush();os.fsync(f.fileno())
def sync_dir(path):
    fd=os.open(path,os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)
def write_new(path,obj):
    with path.open("x") as f: json.dump(obj,f,indent=2,sort_keys=True);f.write("\n");sync(f)
    sync_dir(path.parent)
def reap(proc):
    if proc.poll() is None:
        proc.terminate()
        try:return proc.wait(timeout=5)
        except subprocess.TimeoutExpired: proc.kill();return proc.wait(timeout=5)
    return proc.returncode
def load_plan(path):
    plan=json.loads(path.read_text());validate(plan)
    if plan["frozen_support_sha256"]!=sha(HERE/"plan/frozen_support.py"): raise RuntimeError("frozen source changed after plan validation")
    return plan
def payload(plan,plan_sha):
    frozen=(HERE/"plan/frozen_support.py").read_text();prefix=frozen.split("\ndef main():",1)[0]
    return "SOURCE_SHA="+repr(FROZEN_SHA)+"\n"+prefix+"\n"+(HERE/"remote_support48.py").read_text()+"\nrun_support48(jobs="+json.dumps(plan["jobs"],separators=(",",":"))+", plan_sha256="+repr(plan_sha)+", source_sha256=SOURCE_SHA, workers=15, seconds=14400)\n"
def start(plan_path,run_dir):
    if not Path("/Volumes/data").is_mount() or run_dir.exists(): raise RuntimeError("data unavailable or refusing existing run dir")
    plan=load_plan(plan_path);plan_sha=sha(plan_path);body=payload(plan,plan_sha); body_sha=hashlib.sha256(body.encode()).hexdigest(); permitted={(x["total"],x["m"],x["n"],x["index"],tuple(x["I"])) for x in plan["jobs"]}
    run_dir.mkdir(parents=True);sync_dir(run_dir.parent)
    with (run_dir/"payload.py").open("x") as f:f.write(body);sync(f)
    sync_dir(run_dir)
    command=["ssh",*SSH_OPTIONS,"macbook-math",REMOTE_COMMAND]
    scope="Fresh total-degree 47/48 cap-6 plan only; no total46, no resumed46 coverage, and no general conclusion."
    manifest={"format":"a30-support48-launch-v1","scope":scope,"plan":str(plan_path),"plan_sha256":plan_sha,"frozen_support_sha256":FROZEN_SHA,"remote_sha256":sha(HERE/"remote_support48.py"),"payload_sha256":body_sha,"selected_jobs":len(permitted),"workers":15,"seconds_limit":14400,"ssh_command":command}
    write_new(run_dir/"launch-manifest.json",manifest)
    started=time.time();seen=set();chunks=0;terminal=None;failure=None;rc=None;proc=None
    with (run_dir/"stderr.txt").open("x") as err,(run_dir/"events.jsonl").open("x") as events:
        events.write(json.dumps({"event":"support48_local_started","started_unix":started,"plan_sha256":plan_sha,"payload_sha256":body_sha,"selected_jobs":len(permitted),"scope":scope},separators=(",",":"))+"\n");sync(events)
        try:
            proc=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True,bufsize=1);assert proc.stdin and proc.stdout
            proc.stdin.write(body);proc.stdin.close()
            for line in proc.stdout:
                row=json.loads(line)
                if row.get("event")=="chunk":
                    key=(row["total"],row["m"],row["n"],row["index"],tuple(row["I"]))
                    if key not in permitted or key in seen: raise RuntimeError(f"unplanned/duplicate chunk {key!r}")
                    packed=base64.b64decode(row.pop("gzip_base64"),validate=True)
                    if hashlib.sha256(packed).hexdigest()!=row.get("gzip_sha256"):raise RuntimeError(f"gzip hash {key!r}")
                    raw=gzip.decompress(packed)
                    if hashlib.sha256(raw).hexdigest()!=row.get("raw_sha256"):raise RuntimeError(f"raw hash {key!r}")
                    d=run_dir/f"t{row['total']:02d}"/f"m{row['m']:02d}n{row['n']:02d}";d.mkdir(parents=True,exist_ok=True);sync_dir(d)
                    with (d/f"i{row['index']:06d}.jsonl.gz").open("xb") as f:f.write(packed);sync(f)
                    sync_dir(d);seen.add(key);chunks+=1
                if row.get("event","").startswith("support48_selected_"): terminal=row["event"]
                events.write(json.dumps(row,separators=(",",":"))+"\n");sync(events)
            rc=proc.wait()
        except BaseException as exc:
            failure={"type":type(exc).__name__,"message":str(exc)}
            if proc is not None:rc=reap(proc)
            events.write(json.dumps({"event":"support48_local_error","failure":failure,"chunks_persisted":chunks},separators=(",",":"))+"\n");sync(events)
        finally:
            if proc is not None and proc.poll() is None:rc=reap(proc)
            sync(err)
    complete=rc==0 and terminal=="support48_selected_complete" and seen==permitted
    write_new(run_dir/"exit-receipt.json",{"format":"a30-support48-exit-v1","returncode":rc,"ssh_pid":None if proc is None else proc.pid,"failure":failure,"ended_unix":time.time(),"terminal_event":terminal,"chunks_persisted":chunks,"selected_jobs":len(permitted),"selected_set_complete":complete,"scope":scope})
    return 0 if complete else 1
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--plan",type=Path,required=True);ap.add_argument("--run-dir",type=Path,required=True);a=ap.parse_args();raise SystemExit(start(a.plan,a.run_dir))
if __name__=="__main__":main()
