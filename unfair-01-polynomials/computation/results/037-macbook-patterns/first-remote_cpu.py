"""RAM-only remote discovery. Result ledger travels over stdout to local data disk."""
import base64,gzip,hashlib,json,math,multiprocessing,os,platform,sys,time

def emit(x):
    print(json.dumps(x,separators=(',',':')),flush=True)

def job(pair):
    p,q=pair; begin=time.monotonic(); selected={0}; total=0; rows=[]
    def walk(x):
        if x==q:
            yield tuple(sorted(selected)); return
        if x-p not in selected and x-(q-p) not in selected:
            selected.add(x); yield from walk(x+1); selected.remove(x)
        yield from walk(x+1)
    for T in walk(1):
        total+=1; K=sorted(set(T)|{x+q for x in T}); counts=[0]*p
        for k in K: counts[k%p]+=1
        if 1 not in counts: rows.append(propagate(p,q,K))
    data={'schema':'remote-q-p-propagation-v1','p':p,'q':q,'total_T':total,
          'count':len(rows),'all_unit':sum(not r['unknown_indices'] for r in rows),
          'unknown':sum(bool(r['unknown_indices']) for r in rows),'rows':rows,
          'elapsed_seconds':time.monotonic()-begin,'pid':os.getpid()}
    raw=json.dumps(data,separators=(',',':')).encode(); zipped=gzip.compress(raw,mtime=0)
    return {'event':'result','p':p,'q':q,'count':data['count'],'unknown':data['unknown'],
            'sha256':hashlib.sha256(zipped).hexdigest(),'gzip_base64':base64.b64encode(zipped).decode()}

def run():
    pairs=[(p,q) for q in range(23,27) for p in range(2,q//2+1) if math.gcd(p,q)==1]
    begin=time.monotonic()
    with multiprocessing.get_context('fork').Pool(15) as pool:
        emit({'event':'started','pid':os.getpid(),'worker_pids':[p.pid for p in pool._pool],
              'workers':15,'pairs':pairs,'python':sys.version,'host':platform.node(),
              'source_sha256':SOURCE_SHA,'storage':'RAM only; stdout streamed to Mac mini data'})
        for record in pool.imap_unordered(job,pairs): emit(record)
        pool.close(); pool.join()
    emit({'event':'complete','workers_joined':True,'pairs':len(pairs),'seconds':time.monotonic()-begin})
