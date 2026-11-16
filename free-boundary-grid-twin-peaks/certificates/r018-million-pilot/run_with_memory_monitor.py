#!/usr/bin/env python3
"""Foreground runner with target RSS and global memory-pressure snapshots."""
import json, re, subprocess, sys, threading, time
from pathlib import Path

out=Path(sys.argv[1]); cmd=sys.argv[2:]
proc=subprocess.Popen(cmd); stop=threading.Event()
def sample():
    rss=None; free_pct=None
    try: rss=int(subprocess.check_output(['ps','-o','rss=','-p',str(proc.pid)],text=True).strip())
    except Exception: pass
    try:
        s=subprocess.check_output(['memory_pressure'],text=True,stderr=subprocess.DEVNULL)
        m=re.search(r'memory free percentage:\s*(\d+)%',s,re.I)
        if m: free_pct=int(m.group(1))
    except Exception: pass
    return {'kind':'sample','epoch_seconds':time.time(),'pid':proc.pid,'target_rss_kib':rss,'system_memory_free_pct':free_pct}
with out.open('w') as f:
    f.write(json.dumps({'kind':'metadata','command':cmd,'period_seconds':1,'memory_scope':'global memory_pressure; target ps RSS'})+'\n')
    def loop():
        while not stop.is_set():
            f.write(json.dumps(sample())+'\n');f.flush();stop.wait(1)
    t=threading.Thread(target=loop,daemon=True);t.start();code=proc.wait();stop.set();t.join()
    f.write(json.dumps({'kind':'exit','returncode':code,'epoch_seconds':time.time()})+'\n')
raise SystemExit(code)
