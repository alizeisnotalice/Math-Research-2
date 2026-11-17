"""Rebuild all integer cofactor trees; --mode full also replays every screen.

No success label is accepted as a proof. Trees are rebuilt to their declared
empty layer. Full mode additionally enumerates every Newton mother and checks
all positive-axis and Schur witnesses through the independent implementation.
"""
import argparse,base64,gzip,hashlib,importlib.util,json,math,time,sys,platform
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
from pathlib import Path

def module(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def digest(x):return hashlib.sha256(x).hexdigest()
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()

def screen_init(path,allow):
    global SCREEN,ALLOW
    SCREEN=module(Path(path),'independent_worker');ALLOW=allow
def screen_one(triple):return SCREEN.one(triple,ALLOW)

def main():
    if not __debug__:raise RuntimeError('Certificate verification requires no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--mode',choices=['trees','full'],default='full');ap.add_argument('--workers',type=int,default=1);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();r=a.package;begin=time.monotonic()
    if not 1<=a.workers<=12:raise ValueError('workers must be in 1..12')
    meta=json.loads((r/'manifest.json').read_text())
    for x in meta['files']:
        data=(r/x['path']).read_bytes();assert len(data)==x['bytes'] and digest(data)==x['sha256']
    raw=gzip.decompress((r/'prefix-states.json.gz').read_bytes());assert digest(raw)==meta['prefix_raw_sha256'];prefix=json.loads(raw)
    hgq=module(r/'frozen/verify_s3_hgq.py','hgq');ind=module(r/'frozen/independent_screen.py','independent')
    triples=sorted((T,N,T*T-4*N) for T in range(-16,17) for N in range(-71,72) if T*T-4*N>0 and math.isqrt(T*T-4*N)**2!=T*T-4*N and 144*(T*T-4*N)<(203-12*abs(T))**2)
    assert len(triples)==669
    regenerated=[]
    for T,N,D in triples:
        for c in range(-6,7):
            x=hgq.tree([T,N,D,c,0,0,0],depth=6)
            assert x['status'] in ('unknown_depth','dead_all_prefixes')
            for s in x.get('live_prefixes',[]):
                G,Q=s['G'],s['Q'];regenerated.append([T,N,D,c,G[2]+c,Q[2]+c,*G,*Q])
    assert len(regenerated)==17760 and sorted(regenerated)==sorted(prefix['states'])
    allow={tuple(z[:4]) for z in regenerated};seen=[];count=Counter();mother_count=0;layers=0;depth=0;live=0
    pool=None;futures=None
    if a.mode=='full' and a.workers>1:
        pool=ProcessPoolExecutor(max_workers=a.workers,initializer=screen_init,initargs=(str(r/'frozen/independent_screen.py'),allow))
        futures=[pool.submit(screen_one,tr) for tr in triples]
    parent_source=digest((r/'frozen/parent_screen.py').read_bytes());ind_source=digest((r/'frozen/independent_screen.py').read_bytes())
    for block in meta['blocks']:
        data=(r/block['path']).read_bytes();assert len(data)==block['bytes'] and digest(data)==block['sha256'];n=0
        for line in gzip.decompress(data).splitlines():
            z=json.loads(line);i=z['index'];assert i==len(seen);seen.append(i);n+=1
            docs=[]
            for name in ['parent','independent','certificate']:
                b=z[name+'_utf8'].encode();assert digest(b)==z[name+'_sha256'];docs.append(json.loads(b))
            p,q,cert=docs;assert p['index']==q['index']==cert['index']==i
            assert p['triple']==q['triple']==list(triples[i]) and p['status']==q['status']=='complete'
            assert p['binding']=={'source':parent_source,'prefix':meta['prefix_raw_sha256']}
            assert q['binding']=={'source':ind_source,'prefix':meta['prefix_raw_sha256']}
            assert cert['screen_sha256']==z['parent_sha256']
            packed=base64.b64decode(z['codes_base64'],validate=True);assert digest(packed)==z['codes_gzip_sha256']==p['ledger_sha256'];codes=gzip.decompress(packed)
            assert digest(codes)==p['raw_code_sha256']==q['raw_code_sha256'] and len(codes)==p['counts']['f']
            assert p['counts']==q['counts'] and p['survivors']==q['survivors']
            if a.mode=='full':
                cc,bb,ss=futures[i].result() if futures is not None else ind.one(triples[i],allow);assert cc==p['counts'] and bb==codes and ss==p['survivors']
            count.update(p['counts']);assert [s['row'] for s in cert['records']]==p['survivors'] and len(cert['records'])==cert['mothers']
            for witness in cert['records']:
                k=witness['death'];assert 1<=k<=53
                t=hgq.tree(witness['row'],depth=k,seconds=60)
                assert t['status']=='dead_all_prefixes' and t['depth_reached']==k
                assert digest(canonical(t['trace']))==witness['trace_sha256']
                ml=max(v['live'] for v in t['trace']);assert ml==witness['max_live']
                mother_count+=1;layers+=k;depth=max(depth,k);live=max(live,ml)
        assert n==block['count']
        print(json.dumps({'verified_TN':len(seen),'verified_mothers':mother_count,'mode':a.mode}),flush=True)
    assert seen==list(range(669)) and mother_count==90092 and layers==854628 and depth==53 and live==27
    assert count['f']==27909836 and count['annulus_survived']==90092
    if pool is not None:pool.shutdown(wait=True)
    report={'status':'passed_recomputed_all_90092_empty_trees_and_17760_prefix_states','mode':a.mode,'screen_replayed':a.mode=='full','scope':'Full mode independently replays every necessary screen. Trees mode binds archived screen bytes and rebuilds all cofactor trees; it is not a fresh full screen replay.','TN':669,'mothers':mother_count,'layers':layers,'max_death':depth,'max_live':live,'counts':dict(count),'manifest_sha256':digest((r/'manifest.json').read_bytes()),'verifier_sha256':digest(Path(__file__).read_bytes()),'elapsed_seconds':time.monotonic()-begin}
    report.update({'workers':a.workers,'interpreter':sys.executable,'python_version':platform.python_version(),'argv':sys.argv})
    a.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':main()
