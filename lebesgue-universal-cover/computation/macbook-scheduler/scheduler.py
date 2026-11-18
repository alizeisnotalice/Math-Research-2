"""One account-local Lebesgue scheduler, with a serial GPU lane."""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

POLL = 4
PROBE = 30
START_GRACE = 40
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
ACTIVE = {"starting","running"}

def stamp():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")

def read_json(path):
    return json.loads(path.read_bytes())

def atomic(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp=path.with_name(f".{path.name}.{os.getpid()}")
    with tmp.open("wb") as out:
        out.write(data)
        out.flush()
        os.fsync(out.fileno())
    os.replace(tmp,path)

def save(path,data):
    atomic(path,(json.dumps(data,sort_keys=True,separators=(",",":"))+"\n").encode())

def key_for(task):
    fields=[task[x] for x in ("stage","version","input_sha256","scope")]
    return hashlib.sha256(json.dumps(fields,ensure_ascii=True,separators=(",",":")).encode()).hexdigest()

@contextmanager
def locked(path):
    with path.open("a+") as stream:
        fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        yield

def executable_ps():
    return "/bin/ps" if Path("/bin/ps").exists() else "ps"

def processes():
    text=subprocess.check_output([executable_ps(),"-axo",
        "pid=,ppid=,pgid=,%cpu=,rss=,lstart=,command="],text=True,timeout=8)
    out={}
    for line in text.splitlines():
        p=line.split(None,10)
        if len(p)==11:
            try:
                out[int(p[0])]=dict(ppid=int(p[1]),group=int(p[2]),cpu=float(p[3]),rss=int(p[4]),
                                    start=" ".join(p[5:10]),argv=p[10])
            except ValueError:
                continue
    return out

def descendants(procs, roots):
    ids=set(roots)
    changed=True
    while changed:
        generation={pid for pid,r in procs.items() if r["ppid"] in ids}
        changed=not generation<=ids
        ids.update(generation)
    return ids

def command_identity(procs,pid,needle,start=None):
    proc=procs.get(pid)
    return bool(proc and needle in proc["argv"] and (not start or proc["start"]==start))

def file_result(path, checks):
    if not path.is_file():
        raise ValueError("result file missing")
    minimum=checks.get("minimum_bytes",1)
    if path.stat().st_size<minimum:
        raise ValueError("result file too short")
    kind=checks.get("type","json")
    if kind=="json":
        payload=read_json(path)
        for name,val in checks.get("equals",{}).items():
            probe=payload
            for field in name.split("."):
                probe=probe[field]
            if probe!=val:
                raise ValueError(f"result field {name} does not match")
    elif kind=="jsonl":
        rows=0
        with path.open("rb") as src:
            for line in src:
                if not line.endswith(b"\n"):
                    raise ValueError("incomplete JSONL record")
                row=json.loads(line)
                if not isinstance(row,dict):
                    raise ValueError("JSONL row is not an object")
                if not set(checks.get("row_keys",())) <= set(row):
                    raise ValueError("JSONL row lacks fields")
                rows+=1
        if "rows" in checks and rows!=checks["rows"]:
            raise ValueError(f"expected {checks['rows']} records, found {rows}")
        if rows<checks.get("minimum_rows",1):
            raise ValueError("too few complete JSONL rows")
    elif kind!="file":
        raise ValueError(f"unknown result type {kind}")
    h=hashlib.sha256()
    with path.open("rb") as src:
        for chunk in iter(lambda:src.read(1024*1024),b""):
            h.update(chunk)
    return {"path":str(path),"sha256":h.hexdigest(),"bytes":path.stat().st_size}

def checked_result(run,spec):
    path=(run/spec["path"]).resolve()
    if not path.is_relative_to(run.resolve()):
        raise ValueError("result must be inside its attempt directory")
    return file_result(path,spec)

def render(text,root,run,threads):
    for k,v in {"run":run,"root":root,"threads":threads}.items():
        text=text.replace("{"+k+"}",str(v))
    return text

def worker(root,run):
    with locked(run/"worker.lock"):
        state={"worker_pid":os.getpid(),"started":stamp()}
        save(run/"worker.json",state)
        task=read_json(run/"task.json")
        code=1
        detail=None
        result=None
        try:
            env=os.environ.copy()
            num=task["cpu"]
            env.update({k:str(task.get("math_threads",1)) for k in THREAD_VARS})
            env.update(LEBESGUE_STOP_FILE=str(run/"CHECKPOINT_REQUEST"),
                       LEBESGUE_THREADS=str(num),PYTORCH_ENABLE_MPS_FALLBACK="0")
            command=[render(s,root,run,num) for s in task["argv"]]
            cwd=Path(render(task.get("cwd","{run}"),root,run,num))
            with (run/"stdout.log").open("wb") as out,(run/"stderr.log").open("wb") as err:
                proc=subprocess.Popen(command,cwd=cwd,stdin=subprocess.DEVNULL,stdout=out,stderr=err,env=env)
                save(run/"worker.json",{**state,"child_pid":proc.pid})
                code=proc.wait()
            if code==0:
                spec=task["result"]
                if spec.get("stdout"):
                    atomic(run/spec["path"],(run/"stdout.log").read_bytes())
                result=checked_result(run,spec)
        except (OSError,ValueError,KeyError,TypeError,AttributeError) as exc:
            detail=str(exc)
        save(run/"status.json",{"finished":stamp(),"exit":code,"result":result,"error":detail})
        return 0 if code==0 and result else 1

def sysprobe(argv):
    try:
        return subprocess.check_output(argv,text=True,stderr=subprocess.DEVNULL,timeout=5)
    except (OSError,subprocess.SubprocessError):
        return ""

def measure(previous):
    memory=sysprobe(["/usr/bin/memory_pressure","-Q"])
    swap=sysprobe(["/usr/sbin/sysctl","-n","vm.swapusage"])
    thermal=sysprobe(["/usr/bin/pmset","-g","therm"])
    read=re.search(r"free percentage:\s*(\d+)",memory)
    free=int(read[1]) if read else None
    match=re.search(r"used\s*=\s*([0-9.]+)([MG])",swap)
    used=float(match[1])*(1024 if match[2]=="G" else 1) if match else None
    speed=[int(v) for v in re.findall(r"CPU_Speed_Limit\s*=\s*(\d+)",thermal)]
    increase=used is not None and previous.get("swap_mib") is not None and used>previous["swap_mib"]+32
    growing=previous.get("swap_rises",0)+1 if increase else 0
    critical=(free is not None and free<10) or any(v<50 for v in speed)
    caution=(free is not None and free<22) or growing>=2 or any(v<85 for v in speed)
    level="critical" if critical else "reduced" if caution else "normal"
    calm=previous.get("calm",0)+1 if level=="normal" else 0
    if level=="normal" and previous.get("level","normal")!="normal" and calm<2:
        level="reduced"
    return dict(level=level,memory_free_percent=free,swap_mib=used,swap_rises=growing,
                cpu_speed_limit=min(speed,default=None),calm=calm,measured=stamp())

def validate_task(t):
    for name in ("stage","version","scope"):
        if not isinstance(t.get(name),str) or not t[name].strip():
            raise ValueError(f"missing {name}")
    if not re.fullmatch(r"[0-9a-f]{64}",t.get("input_sha256","")):
        raise ValueError("input_sha256 must be 64 lowercase hex characters")
    if t.get("kind") not in ("cpu","gpu") or type(t.get("cpu")) is not int or not 1<=t["cpu"]<=128:
        raise ValueError("specify kind and positive CPU reservation")
    if type(t.get("priority",0)) is not int:
        raise ValueError("priority must be an integer")
    spec=t["result"]
    result=Path(spec["path"])
    if result.is_absolute() or ".." in result.parts or result.name in ("status.json","worker.json","task.json"):
        raise ValueError("invalid result path")
    if not all(re.fullmatch(r"[0-9a-f]{64}",s) for s in t.get("after",())):
        raise ValueError("dependencies must be task keys")
    if t.get("driver")=="seven":
        if t["kind"]!="gpu" or type(t["epochs"]) is not int or not 1<=t["epochs"]<=50000:
            raise ValueError("invalid native GPU task or epoch count")
        if type(t["seed"]) is not int or not 0<=t["seed"]<2**32:
            raise ValueError("seed outside UInt32")
        if spec["path"]!="result.json" or not Path(t["input"]).is_absolute():
            raise ValueError("native GPU input must be absolute and output is result.json")
    elif not t["argv"] or not all(isinstance(x,str) for x in t["argv"]) or not t["argv"][0].startswith("/"):
        raise ValueError("argv must start with an absolute executable")
    return key_for(t)

class Scheduler:
    def __init__(self,root):
        self.root=root
        self.ledger=root/"ledger.json"
        self.book=read_json(self.ledger) if self.ledger.exists() else {"schema":1,"tasks":{},"external":{}}
        self.tasks=self.book["tasks"]
        self.external=self.book["external"]
        self.seen=set()
        self.issues={}
        self.probe_time=0
        self.pressure={"level":"normal"}
        self.stop=False
        self.last_emit=0
        self.gpu_child=None
        self.children=[]

    def write(self):
        save(self.ledger,self.book)

    def ingest(self):
        for path in sorted((self.root/"inbox").glob("*.json")):
            raw=path.read_bytes()
            ticket=(path.name,hashlib.sha256(raw).hexdigest())
            if ticket in self.seen:
                continue
            try:
                rows=json.loads(raw)["tasks"]
                keys=[validate_task(t) for t in rows]
                if len(set(keys))!=len(keys):
                    raise ValueError("duplicate task key inside one manifest")
                if any(k in self.tasks and self.tasks[k]["task"]!=t for k,t in zip(keys,rows,strict=True)):
                    raise ValueError("a task key has a different definition")
                for k,t in zip(keys,rows,strict=True):
                    existing=self.tasks.get(k)
                    if existing is None:
                        self.tasks[k]={"task":t,"state":"queued","queued":stamp(),"manifest":path.name}
                for k in json.loads(raw).get("retry",[]):
                    item=self.tasks.get(k)
                    if item is not None and item["state"] in ("failed","interrupted") and not item.get("last_retry")==ticket[1]:
                        item["last_retry"]=ticket[1]
                        item["history"]=item.get("history",[])+[{a:item.get(a) for a in ("state","run","finished","error","exit")}]
                        item.update(state="queued",queued=stamp())
                self.issues.pop(path.name,None)
            except (ValueError,KeyError,TypeError) as exc:
                self.issues[path.name]=str(exc)
            self.seen.add(ticket)
        live=self.root/"external.json"
        if live.exists():
            for t in read_json(live)["jobs"]:
                name=t["id"]
                if name not in self.external:
                    self.external[name]={"job":t,"state":"new","descendant_starts":{}}
                elif self.external[name]["job"]!=t:
                    self.issues[name]="external identifier already bound"

    def external_status(self,procs):
        roots=set()
        gpu=[]
        cpus=0
        for name,item in self.external.items():
            spec=item["job"]
            pid=spec["pid"]
            active=command_identity(procs,pid,spec["command_contains"],item.get("started") or spec.get("start"))
            if active:
                item["started"]=procs[pid]["start"]
                fresh=descendants(procs,[pid])
                item["descendant_starts"].update({str(k):procs[k]["start"] for k in fresh if k!=pid and k in procs})
            # Remember already observed children if their parent exits first.
            remaining={k:s for k,s in item["descendant_starts"].items()
                       if command_identity(procs,int(k),"",s)}
            item["descendant_starts"]=remaining
            if active or remaining:
                item["state"]="attached"
                roots.update(map(int,remaining))
                if active:
                    roots.add(pid)
                cpus+=spec["cpu"]
                if spec.get("gpu"):
                    gpu.append(name)
                continue
            if item["state"] not in ("complete","ended_unverified","failed"):
                item["ended"]=stamp()
                result=spec.get("result")
                try:
                    if not result:
                        raise ValueError("legacy job had no completion contract")
                    if "exit_path" in spec and Path(spec["exit_path"]).read_text().strip()!="0":
                        item["state"]="failed"
                        continue
                    item["receipt"]=file_result(Path(result["path"]),result)
                    item["state"]="complete"
                except (OSError,KeyError,TypeError,ValueError) as exc:
                    item.update(state="ended_unverified",error=str(exc))
        return cpus,gpu,roots

    def native_pid(self,procs):
        info=self.root/"gpu-spool"/"service.pid"
        if not info.is_file():
            return None
        try:
            pid=int(info.read_text().strip())
        except ValueError:
            return None
        needle=str(self.root/"gpu-spool")
        return pid if command_identity(procs,pid,needle) else None

    def reconcile(self,procs):
        owned=set()
        for k,row in self.tasks.items():
            if row["state"] not in ACTIVE:
                continue
            run=Path(row["run"])
            status=run/"status.json"
            if status.exists():
                recap=read_json(status)
                try:
                    if recap["exit"]==0:
                        result=checked_result(run,row["task"]["result"])
                        if recap.get("result") and result!=recap["result"]:
                            raise ValueError("result changed after worker verification")
                        row.update(state="complete",result=result,finished=recap["finished"])
                    else:
                        row.update(state="failed",error=recap.get("error"),exit=recap["exit"],
                                   finished=recap["finished"])
                except (OSError,ValueError,KeyError,TypeError) as exc:
                    row.update(state="failed",error=str(exc),finished=stamp())
                continue
            native=row["task"].get("driver")=="seven"
            if native:
                pid=self.native_pid(procs)
            else:
                info=run/"worker.json"
                pid=read_json(info)["worker_pid"] if info.exists() else row.get("pid")
                if pid and not command_identity(procs,pid,str(run)):
                    pid=None
            if pid:
                row["pid"]=pid
                row["state"]="running"
                owned.add(pid)
            elif time.time()-row["claimed_at"]>START_GRACE:
                row.update(state="interrupted",error="process ended without a complete receipt",finished=stamp())
        return owned

    def claim(self,k,row,procs):
        run=self.root/"runs"/k/f"attempt-{1+len(row.get('history',[])):03d}"
        run.mkdir(parents=True)
        task=row["task"]
        save(run/"task.json",task)
        row.update(state="starting",run=str(run),claimed=stamp(),claimed_at=time.time())
        self.write()
        try:
            if task.get("driver")=="seven":
                spool=self.root/"gpu-spool"
                cfg=read_json(self.root/"config.json")
                inp=Path(task["input"])
                if hashlib.sha256(inp.read_bytes()).hexdigest()!=task["input_sha256"]:
                    raise ValueError("native GPU input SHA-256 changed")
                if self.native_pid(procs) is None:
                    if self.gpu_child is None or self.gpu_child.poll() is not None:
                        with (spool/"service.log").open("ab") as log:
                            self.gpu_child=subprocess.Popen([str(self.root/"resident-gpu"),
                                cfg["metal_source"],str(spool)],stdin=subprocess.DEVNULL,
                                stdout=log,stderr=log,start_new_session=True)
                    row["pid"]=self.gpu_child.pid
                else:
                    row["pid"]=self.native_pid(procs)
                fields=[json.dumps(str(inp)),json.dumps(str(run)),str(task["epochs"]),str(task["seed"])]
                atomic(spool/(k+".job"),(" ".join(fields)+"\n").encode())
            else:
                proc=subprocess.Popen([sys.executable,"-B",str(Path(__file__).resolve()),"worker",str(self.root),str(run)],
                    stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True)
                row["pid"]=proc.pid
                self.children.append(proc)
            self.write()
        except (OSError,ValueError,KeyError) as exc:
            row.update(state="failed",error=str(exc),finished=stamp())
            self.write()

    def step(self):
        self.children=[p for p in self.children if p.poll() is None]
        self.ingest()
        procs=processes()
        own=self.reconcile(procs)
        external_cpus, external_gpu, observed=self.external_status(procs)
        owned_tree=descendants(procs,{os.getpid(),*own,*observed})
        foreign=sum(row["cpu"] for pid,row in procs.items() if pid not in owned_tree)
        cfg=read_json(self.root/"config.json")
        if time.monotonic()-self.probe_time>=PROBE:
            self.pressure=measure(self.pressure)
            self.probe_time=time.monotonic()
            save(self.root/"gpu-spool"/"pressure.json",self.pressure)
        maximum=max(0,cfg["cpu_budget"]-math.ceil(max(0,foreign-20)/100))
        if self.pressure["level"]=="reduced":
            maximum//=2
        elif self.pressure["level"]=="critical":
            maximum=0
        active=[row for row in self.tasks.values() if row["state"] in ACTIVE]
        gpu_busy=bool(external_gpu or any(t["task"]["kind"]=="gpu" for t in active))
        free=max(0,maximum-external_cpus-sum(r["task"]["cpu"] for r in active))
        stopping=self.stop or (self.root/"STOP").exists()
        if stopping or self.pressure["level"]=="critical":
            eligible=[r for r in active if r["task"].get("checkpoint")]
            if not stopping:
                eligible=sorted(eligible,key=lambda r:sum(procs[p]["rss"] for p in descendants(procs,[r.get("pid")]) if p in procs),
                                reverse=True)[:1]
            for r in eligible:
                mark=Path(r["run"])/"CHECKPOINT_REQUEST"
                if not mark.exists():
                    atomic(mark,(stamp()+"\n").encode())
        if not stopping:
            queue=sorted(((k,t) for k,t in self.tasks.items() if t["state"]=="queued"),
                key=lambda x:(x[1]["task"]["kind"]!="cpu",-x[1]["task"].get("priority",0),x[1]["queued"],x[0]))
            ready=[(k,r) for k,r in queue if all(self.tasks.get(key,{}).get("state")=="complete" for key in r["task"].get("after",()))]
            reserve=next((r["task"]["cpu"] for _,r in ready if r["task"]["kind"]=="gpu" and r["task"]["cpu"]<=free),0) if not gpu_busy else 0
            for k,row in ready:
                t=row["task"]
                if t["cpu"]>free-(reserve if t["kind"]=="cpu" else 0):
                    continue
                if t["kind"]=="gpu" and gpu_busy:
                    continue
                self.claim(k,row,procs)
                if row["state"] in ACTIVE:
                    free-=t["cpu"]
                    gpu_busy|=t["kind"]=="gpu"
        self.write()
        todo={resource:sum(r["state"]=="queued" and r["task"]["kind"]==resource for r in self.tasks.values())
              for resource in ("cpu","gpu")}
        stats={state:sum(r["state"]==state for r in self.tasks.values())
               for state in ("queued","starting","running","complete","failed","interrupted")}
        save(self.root/"heartbeat.json",dict(pid=os.getpid(),host=os.uname().nodename,time=stamp(),
             pressure=self.pressure,logical_cpu=os.cpu_count(),configured_cpu_budget=cfg["cpu_budget"],effective_budget=maximum,
             external_cpu=external_cpus,foreign_cpu_percent=round(foreign,1),available_cpu=free,
             external_gpu=external_gpu,native_gpu_pid=self.native_pid(procs),queue=todo,counts=stats,
             refill={"cpu":todo["cpu"]<=maximum,"gpu":todo["gpu"]<=2},stopping=stopping,manifest_errors=self.issues,
             running=[{"key":k,"pid":r.get("pid"),"scope":r["task"]["scope"],"kind":r["task"]["kind"]}
                      for k,r in self.tasks.items() if r["state"] in ACTIVE]))
        if stopping and not any(r["state"] in ACTIVE for r in self.tasks.values()):
            atomic(self.root/"gpu-spool"/"STOP",b"drain\n")
            return False
        return True

    def serve(self):
        def stopped(*_):
            self.stop=True
        signal.signal(signal.SIGTERM,stopped)
        signal.signal(signal.SIGINT,stopped)
        while True:
            try:
                if not self.step():
                    break
            except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as exc:
                print(stamp(),type(exc).__name__,str(exc),flush=True)
            time.sleep(POLL)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=("serve","detach","worker","key"))
    p.add_argument("root",type=Path)
    p.add_argument("run",type=Path,nargs="?")
    a=p.parse_args()
    root=a.root.resolve()
    if a.action=="key":
        print(validate_task(read_json(root)))
        return 0
    for name in ("inbox","runs","inputs","gpu-spool"):
        (root/name).mkdir(parents=True,exist_ok=True)
    if a.action=="worker":
        return worker(root,a.run)
    if a.action=="detach":
        if (root/"STOP").exists():
            raise SystemExit("Project STOP is present")
        with (root/"service.log").open("ab") as log:
            child=subprocess.Popen([sys.executable,"-B",str(Path(__file__).resolve()),"serve",str(root)],
                 stdin=subprocess.DEVNULL,stdout=log,stderr=log,start_new_session=True)
        print(child.pid)
        return 0
    try:
        with locked(root/"scheduler.lock"):
            (root/"gpu-spool"/"STOP").unlink(missing_ok=True)
            atomic(root/"scheduler.pid",f"{os.getpid()}\n".encode())
            Scheduler(root).serve()
        return 0
    except BlockingIOError:
        print("Lebesgue scheduler already holds the lease",file=sys.stderr)
        return 73

if __name__=="__main__":
    raise SystemExit(main())
