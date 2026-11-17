#!/usr/bin/env python3
"""Remote-only controller for the GPU52 continuation.

It never writes task data on the remote host: all JSON records pass through
stdout to the local bridge.  It accepts a task only after validating its exact
resume boundary and it always reaps any started WolframKernel.
"""
import base64
import gc
import hashlib
import json
import os
import subprocess
import sys
import time

all_data = json.loads(base64.b64decode(INPUT_B64))
search = base64.b64decode(SEARCH_B64).decode()
SECONDS_LIMIT = int(os.environ.get("GPU52_SECONDS_LIMIT", "14400"))
SMOKE = os.environ.get("GPU52_SMOKE", "0") == "1"


def emit(row):
    print(json.dumps(row, separators=(",", ":")), flush=True)


started = time.monotonic()
covered = 0
launched = 0
skipped = 0
emit({"event": "controller_started", "pid": os.getpid(), "python": sys.version,
      "tasks": len(all_data["tasks"]), "seconds_limit": SECONDS_LIMIT,
      "smoke": SMOKE, "storage": "RAM only; no remote task files"})
for task in all_data["tasks"]:
    task_id, total = int(task["id"]), int(task["masks"])
    start = task.get("resume_start_mask")
    if not isinstance(start, int) or not 0 <= start <= total:
        emit({"event": "invalid_resume_start", "task": task_id, "start": start, "masks": total})
        break
    if start == total:
        skipped += 1
        covered += 1
        emit({"event": "task_skipped_completed", "task": task_id, "start": start,
              "masks": total, "reason": "old_stream_audited_prefix"})
        continue
    remaining = SECONDS_LIMIT - (time.monotonic() - started)
    if remaining <= 1:
        emit({"event": "controller_deadline_before_task", "task": task_id})
        break
    run_task = dict(task)
    if SMOKE:
        # Explicitly bounded test.  Its record must never be merged into the
        # production coverage ledger, which still begins at resume_start_mask.
        run_task["masks"] = min(total, start + int(all_data["batch"]))
    data = {**all_data, "tasks": [run_task]}
    wl = "data=ImportString[" + json.dumps(json.dumps(data, separators=(",", ":"))) + ',"RawJSON"];\n' + search
    if len(wl.encode()) >= 245000:
        emit({"event": "program_too_large", "task": task_id, "bytes": len(wl.encode())})
        break
    p = subprocess.Popen(["/Applications/Wolfram.app/Contents/MacOS/WolframKernel", "-noprompt", "-run", wl], stdin=subprocess.DEVNULL)
    launched += 1
    emit({"event": "kernel_started", "task": task_id, "pid": p.pid,
          "resume_first_mask": start, "planned_end": int(run_task["masks"]),
          "program_sha256": hashlib.sha256(wl.encode()).hexdigest(), "program_bytes": len(wl.encode())})
    timed_out = False
    try:
        rc = p.wait(timeout=remaining)
    except subprocess.TimeoutExpired:
        timed_out = True
        p.terminate()
        try:
            rc = p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            p.kill()
            rc = p.wait()
    finally:
        if p.poll() is None:
            p.kill()
            rc = p.wait()
    emit({"event": "kernel_reaped", "task": task_id, "pid": p.pid, "returncode": rc, "timed_out": timed_out})
    del data, wl
    gc.collect()
    if rc != 0 or timed_out:
        break
    # A launcher return code is not a completion certificate; the local
    # verifier must demand full contiguous batch records and terminal events.
    covered += 1
    if SMOKE:
        emit({"event": "controller_smoke_complete", "task": task_id, "all_started_kernels_reaped": True})
        break
else:
    pass
terminal = "controller_complete" if covered == len(all_data["tasks"]) and not SMOKE else "controller_partial"
emit({"event": terminal, "covered_task_rows": covered, "tasks": len(all_data["tasks"]),
      "launched": launched, "skipped": skipped, "seconds": time.monotonic()-started,
      "all_started_kernels_reaped": True, "smoke": SMOKE})
