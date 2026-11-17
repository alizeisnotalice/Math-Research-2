"""Injected below the frozen classifier; use selected missing total46 jobs only."""
import base64, io, multiprocessing, platform


def compute_selected(job):
    m, n, index, I = job["m"], job["n"], job["index"], tuple(job["I"])
    assert m + n == 46 and 1 <= m <= 23 and I[0] == 0 and I[-1] == m
    counts = dict(raw=0, canonical=0, excluded=0, fair=0, unknown=0)
    raw_stream = io.BytesIO(); digest = hashlib.sha256(); unknown = []
    with gzip.GzipFile(fileobj=raw_stream, mode="wb", compresslevel=1, mtime=0) as z:
        for J in supports(n, 6):
            counts["raw"] += 1
            if not canonical(I, J):
                continue
            status, witness = classify(I, J)
            line = (json.dumps([I, J, status, witness], separators=(",", ":")) + "\n").encode()
            z.write(line); digest.update(line); counts["canonical"] += 1; counts[status] += 1
            if status == "unknown":
                unknown.append([I, J])
    packed = raw_stream.getvalue()
    return {"event": "chunk", "m": m, "n": n, "index": index, "I": I, "counts": counts,
            "unknown": unknown, "raw_sha256": digest.hexdigest(), "gzip_sha256": hashlib.sha256(packed).hexdigest(),
            "gzip_base64": base64.b64encode(packed).decode()}


def emit_selected(row):
    print(json.dumps(row, separators=(",", ":")), flush=True)


def run_selected(jobs, plan_sha256, source_sha256, workers=15, seconds=14400):
    assert workers == 15 and seconds == 14400
    begin = time.monotonic(); pending = []; cursor = 0; timed_out = False; completed = 0
    emit_selected({"event": "resume46_remote_started", "pid": os.getpid(), "worker_pids": [], "workers": workers,
                   "python": sys.version, "host": platform.node(), "frozen_source_sha256": source_sha256,
                   "plan_sha256": plan_sha256, "selected_jobs": len(jobs), "total": 46, "cap": 6,
                   "seconds_limit": seconds, "storage": "RAM only; chunks stream to local data"})
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        emit_selected({"event": "resume46_remote_workers", "worker_pids": [p.pid for p in pool._pool]})
        while True:
            while cursor < len(jobs) and len(pending) < 2 * workers and time.monotonic() - begin < seconds:
                pending.append(pool.apply_async(compute_selected, (jobs[cursor],))); cursor += 1
            if cursor < len(jobs) and time.monotonic() - begin >= seconds:
                timed_out = True
            ready = [x for x in pending if x.ready()]
            for item in ready:
                emit_selected(item.get()); pending.remove(item); completed += 1
            if not pending and (cursor == len(jobs) or timed_out):
                break
            if not ready:
                time.sleep(.01)
        pool.close(); pool.join()
    status = "resume46_selected_complete" if cursor == len(jobs) and not timed_out else "resume46_selected_partial_timeout"
    emit_selected({"event": status, "workers_joined": True, "submitted": cursor, "finished": completed,
                   "selected_jobs": len(jobs), "complete_selected_set": status == "resume46_selected_complete",
                   "seconds": time.monotonic() - begin,
                   "scope": "Only missing total-degree-46 chunks; does not certify the interrupted original 45/46 batch."})
