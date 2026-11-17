#!/usr/bin/env python3
"""Bounded parallel executor for round014 exact capacity-DP tasks."""
import argparse, concurrent.futures, hashlib, json, os, queue, re, signal, subprocess, threading, time
from pathlib import Path

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def idle():
    try:
        x=subprocess.check_output(['top','-l','1','-n','0'],text=True,stderr=subprocess.DEVNULL,timeout=5)
        m=re.search(r'([0-9.]+)% idle',x)
        return float(m.group(1)) if m else None
    except Exception:return None

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--tasks',type=Path,required=True);p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--results',type=Path,required=True);p.add_argument('--summary',type=Path,required=True)
    p.add_argument('--workers',type=int,default=15);p.add_argument('--per-task-seconds',type=float,default=90)
    p.add_argument('--batch-seconds',type=float,default=600);p.add_argument('--sample-seconds',type=float,default=2)
    a=p.parse_args(); a.binary=a.binary.resolve(); tasks=[json.loads(x) for x in a.tasks.read_text().splitlines() if x]
    began=time.monotonic();deadline=began+a.batch_seconds;running={};records=[];samples=[];lock=threading.Lock();stop=threading.Event()
    def monitor():
        while not stop.is_set():
            with lock:pids=list(running)
            child=0.0
            if pids:
                try: child=sum(float(x) for x in subprocess.check_output(['ps','-p',','.join(map(str,pids)),'-o','%cpu='],text=True,stderr=subprocess.DEVNULL).split())
                except Exception:pass
            samples.append({'seconds':time.monotonic()-began,'cpu_idle_percent':idle(),'sum_child_cpu_percent':child,'active_children':len(pids)})
            stop.wait(a.sample_seconds)
    def run(task):
        if time.monotonic()>=deadline:return {'id':task['id'],'status':'not_started_batch_deadline'}
        cmd=[str(a.binary),str(task['lift_dimension'])]+[str(x) for x in task['lift_basis']];t=time.monotonic()
        q=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)
        with lock:running[q.pid]=task['id']
        try:out,err=q.communicate(timeout=min(a.per_task_seconds,max(.1,deadline-time.monotonic())))
        except subprocess.TimeoutExpired:
            os.killpg(q.pid,signal.SIGTERM);out,err=q.communicate()
            with lock:running.pop(q.pid,None)
            return {'id':task['id'],'status':'timeout','elapsed_seconds':time.monotonic()-t,'task':task,'stdout':out,'stderr':err}
        with lock:running.pop(q.pid,None)
        r={'id':task['id'],'elapsed_seconds':time.monotonic()-t,'task':task,'returncode':q.returncode,'stdout':out,'stderr':err}
        if q.returncode:r['status']='nonzero_exit';return r
        try:
            z=json.loads(out);r.update({'status':'ok','D':z['D'],'kappa':z['kappa'],'partition':z['partition'],'upper':task['upper'],'violates_upper':z['kappa']>task['upper']})
        except Exception as e:r.update({'status':'malformed_json','parse_error':repr(e)})
        return r
    pending=queue.Queue()
    for task in tasks: pending.put(task)
    a.results.parent.mkdir(parents=True,exist_ok=True);m=threading.Thread(target=monitor,daemon=True);m.start()
    with open(a.results,'w') as f:
        def consume():
            while time.monotonic()<deadline:
                try: task=pending.get_nowait()
                except queue.Empty: return
                record=run(task)
                with lock:
                    records.append(record);f.write(json.dumps(record,separators=(',',':'))+'\n');f.flush()
        with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as ex:
            fs=[ex.submit(consume) for _ in range(a.workers)]
            for future in fs: future.result()
    stop.set();m.join(timeout=10);ok=[x for x in records if x['status']=='ok'];bad=[x for x in records if x['status']!='ok'];ce=[x for x in ok if x['violates_upper']]
    summary={'status':'complete','tasks_total':len(tasks),'records':len(records),'unstarted_tasks':pending.qsize(),'successful':len(ok),'failures':bad,'counterexamples':ce,
      'boundary_records':sorted(ok,key=lambda x:(x['kappa']-x['upper'],x['kappa']),reverse=True)[:20],'workers':a.workers,'per_task_seconds':a.per_task_seconds,'batch_seconds':a.batch_seconds,'elapsed_seconds':time.monotonic()-began,'cpu_samples':samples,'tasks_sha256':sha256(a.tasks),'binary_sha256':sha256(a.binary),'results_sha256':sha256(a.results)}
    a.summary.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({'status':'complete','successful':len(ok),'failures':len(bad),'counterexamples':len(ce),'elapsed_seconds':summary['elapsed_seconds']}))
if __name__=='__main__':main()
