"""Independent mask coverage and replay; run entirely in remote RAM."""
import base64,gzip,hashlib,json,multiprocessing,os,sys,time

def verify_job(item):
    global Q
    name,encoded=item; raw=base64.b64decode(encoded); d=json.loads(gzip.decompress(raw))
    p,Q=d['p'],d['q']; started=time.monotonic(); expected=set(); eligible=0
    for mask in range(1,1<<Q,2):
        g=mask|(mask<<Q)
        if g&(g>>p): continue
        eligible+=1
        if no_singleton(mask,p):expected.add(k_from_mask(mask))
    seen=set(); events=Counter(); unknown=0
    for row in d['rows']:
        K=tuple(row['G_support'])
        if row['p']!=p or row['q']!=Q or K not in expected or K in seen:raise ValueError('identity/coverage')
        seen.add(K); events.update(replay(row,(p,K)));unknown+=bool(row['unknown_indices'])
    if seen!=expected or eligible!=d['total_T'] or len(seen)!=d['count'] or unknown!=d['unknown'] or len(seen)-unknown!=d['all_unit']:raise ValueError('complete coverage/count mismatch')
    return {'file':name,'sha256':hashlib.sha256(raw).hexdigest(),'p':p,'q':Q,'masks':1<<(Q-1),'eligible':eligible,'count':len(seen),'unknown':unknown,'events':dict(events),'seconds':time.monotonic()-started,'pid':os.getpid(),'verified':True}

def run_verification():
    begin=time.monotonic()
    with multiprocessing.get_context('fork').Pool(15) as pool:
        print(json.dumps({'event':'started','pid':os.getpid(),'workers':[p.pid for p in pool._pool],'python':sys.version}),flush=True)
        for r in pool.imap_unordered(verify_job,INPUTS): print(json.dumps(r),flush=True)
        pool.close();pool.join()
    print(json.dumps({'event':'complete','workers_joined':True,'seconds':time.monotonic()-begin}),flush=True)
