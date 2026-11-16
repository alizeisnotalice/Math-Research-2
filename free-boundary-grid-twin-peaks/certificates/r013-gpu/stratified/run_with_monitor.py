#!/usr/bin/env python3
"""Run one foreground command while sampling whole-machine AGX counters and RSS.

The AGX values come from IOKit's global PerformanceStatistics; they include
other clients and are not per-process utilization.  powermetrics is not used,
because it requires superuser privilege on this host.
"""
import json
from pathlib import Path
import re
import subprocess
import sys
import threading
import time

out = Path(sys.argv[1])
cmd = sys.argv[2:]
if not cmd:
    raise SystemExit("usage: run_with_monitor.py telemetry.jsonl command ...")

proc = subprocess.Popen(cmd)
stop = threading.Event()

def sample():
    rss = None
    try:
        rss = int(subprocess.check_output(["ps", "-o", "rss=", "-p", str(proc.pid)], text=True).strip())
    except (subprocess.CalledProcessError, ValueError):
        pass
    agx = {"device_utilization_pct": None, "renderer_utilization_pct": None, "tiler_utilization_pct": None}
    try:
        text = subprocess.check_output(["ioreg", "-r", "-d1", "-c", "AGXAccelerator"], text=True, stderr=subprocess.DEVNULL)
        for key, name in [("device_utilization_pct", "Device Utilization %"), ("renderer_utilization_pct", "Renderer Utilization %"), ("tiler_utilization_pct", "Tiler Utilization %")]:
            m = re.search(r'"' + re.escape(name) + r'"=(\d+)', text)
            if m: agx[key] = int(m.group(1))
    except subprocess.CalledProcessError:
        pass
    return {"epoch_seconds": time.time(), "pid": proc.pid, "rss_kib": rss, **agx}

with out.open("w") as f:
    f.write(json.dumps({"kind": "metadata", "command": cmd, "period_seconds": 1,
                        "agx_scope": "whole-machine IOKit PerformanceStatistics, nonexclusive",
                        "rss_scope": "target process ps rss", "powermetrics": "unavailable without superuser"}) + "\n")
    def loop():
        while not stop.is_set():
            f.write(json.dumps({"kind": "sample", **sample()}) + "\n"); f.flush()
            stop.wait(1)
    t = threading.Thread(target=loop, daemon=True); t.start()
    code = proc.wait()
    stop.set(); t.join()
    f.write(json.dumps({"kind": "exit", "returncode": code, "epoch_seconds": time.time()}) + "\n")
raise SystemExit(code)
