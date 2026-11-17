"""Injected after the frozen classifier; fresh total-47/48 chunks only."""
import base64, io, multiprocessing, platform

def compute_support48(job):
    total,m,n,index,I=job["total"],job["m"],job["n"],job["index"],tuple(job["I"])
    # The local bridge has independently validated the complete I/index table.
    # Keep this worker guard O(1); rebuilding the entire table for every one of
    # 100005 chunks would be a quadratic accidental cost.
    assert total in (47,48) and m+n==total and 1<=m<=total//2 and index>=0 and I[0]==0 and I[-1]==m
    counts=dict(raw=0,canonical=0,excluded=0,fair=0,unknown=0); buf=io.BytesIO(); digest=hashlib.sha256(); unknown=[]
    with gzip.GzipFile(fileobj=buf,mode="wb",compresslevel=1,mtime=0) as z:
        for J in supports(n,6):
            counts["raw"]+=1
            if not canonical(I,J): continue
            status,witness=classify(I,J); line=(json.dumps([I,J,status,witness],separators=(",",":"))+"\n").encode()
            z.write(line); digest.update(line); counts["canonical"]+=1; counts[status]+=1
            if status=="unknown": unknown.append([I,J])
    packed=buf.getvalue()
    return {"event":"chunk","total":total,"m":m,"n":n,"index":index,"I":I,"counts":counts,"unknown":unknown,"raw_sha256":digest.hexdigest(),"gzip_sha256":hashlib.sha256(packed).hexdigest(),"gzip_base64":base64.b64encode(packed).decode()}

def emit48(row): print(json.dumps(row,separators=(",",":")),flush=True)

def run_support48(jobs,plan_sha256,source_sha256,workers=15,seconds=14400):
    assert workers==15 and seconds==14400
    permitted={(x["total"],x["m"],x["n"],x["index"],tuple(x["I"])) for x in jobs}; assert len(permitted)==len(jobs)
    begin=time.monotonic(); pending=[]; cursor=finished=0; timed_out=False
    emit48({"event":"support48_remote_started","pid":os.getpid(),"worker_pids":[],"workers":workers,"python":sys.version,"host":platform.node(),"frozen_source_sha256":source_sha256,"plan_sha256":plan_sha256,"selected_jobs":len(jobs),"total_min":47,"total_max":48,"cap":6,"seconds_limit":seconds,"storage":"worker RAM only; completed chunks stream to local data ledger"})
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        emit48({"event":"support48_remote_workers","worker_pids":[p.pid for p in pool._pool]})
        while True:
            while cursor<len(jobs) and len(pending)<2*workers and time.monotonic()-begin<seconds:
                pending.append(pool.apply_async(compute_support48,(jobs[cursor],))); cursor+=1
            if cursor<len(jobs) and time.monotonic()-begin>=seconds: timed_out=True
            ready=[x for x in pending if x.ready()]
            for item in ready: emit48(item.get());pending.remove(item);finished+=1
            if not pending and (cursor==len(jobs) or timed_out): break
            if not ready: time.sleep(.01)
        pool.close();pool.join()
    event="support48_selected_complete" if cursor==len(jobs) and not timed_out else "support48_selected_partial_timeout"
    emit48({"event":event,"workers_joined":True,"submitted":cursor,"finished":finished,"selected_jobs":len(jobs),"complete_selected_set":event=="support48_selected_complete","seconds":time.monotonic()-begin,"scope":"Fresh total-degree 47/48 selected plan only; no total46 or general conclusion."})
