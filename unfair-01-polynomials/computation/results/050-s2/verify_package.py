"""Portable exact S2 certificate replay; full is the default mode."""
import argparse,base64,contextlib,gzip,hashlib,importlib.util,json,platform,sys,time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
def digest(b):return hashlib.sha256(b).hexdigest()
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def module(p,name):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def screen_init(path):
    global SCREEN
    SCREEN=module(Path(path),'independent_worker')
def screen_one(core):return SCREEN.screen_one(core)
def counts_normalized(c,codes):
    hist=Counter(codes)
    for k in range(4,8):assert c.get(f'reject{k}',0)==hist[k]
    return Counter({k.replace('axis_','axis'):v for k,v in c.items() if not k.startswith('reject')})
def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--mode',choices=['full','trees'],default='full');ap.add_argument('--workers',type=int,default=1);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();r=args.package;begin=time.monotonic();assert 1<=args.workers<=12
    meta=json.loads((r/'manifest.json').read_text())
    for x in meta['files']:
        b=(r/x['path']).read_bytes();assert len(b)==x['bytes'] and digest(b)==x['sha256']
    raw=gzip.decompress((r/'prefix-ledger.json.gz').read_bytes());assert digest(raw)==meta['prefix_raw_sha256'];prefix=json.loads(raw)
    pre=module(r/'frozen/prefix_independent.py','prefix_independent');hgq=module(r/'frozen/verify_s2_hgq.py','hgq')
    core=pre.corrected_core();expected_core={(z['T'],z['N'],z['Delta'],z['b']) for z in prefix['core_records']}
    assert len(core)==len(prefix['core_records'])==4913 and core==expected_core
    low=pre.independently_surviving(core);fields=('T','N','Delta','b','c','A1','B1','G1','G2','Q1','Q2','c3');expected_low={tuple(z[k] for k in fields) for z in prefix['survivors']}
    assert len(low)==len(prefix['survivors'])==395 and low==expected_low
    cores=sorted({tuple(z[:5]) for z in low});assert len(cores)==359
    regenerated=[]
    for c in cores:
        tr=hgq.tree([*c,0,0,0],depth=3)
        regenerated.extend([*c,*z['G'],*z['Q']] for z in tr.get('live_prefixes',[]))
    want=[[*[z[k] for k in fields[:5]],1,z['G1'],z['G2'],0,1,z['Q1'],z['Q2'],0] for z in prefix['survivors']]
    assert sorted(regenerated)==sorted(want);print(json.dumps({'verified_core_box':4913,'verified_prefixes':395,'seconds':time.monotonic()-begin}),flush=True)
    parent_source=digest((r/'frozen/parent_screen.py').read_bytes());ind_source=digest((r/'frozen/independent_screen.py').read_bytes());seen=[];mother_seen=set();total=Counter();layers=0;maxdeath=0;maxlive=0;depths=Counter()
    poolctx=ProcessPoolExecutor(max_workers=args.workers,initializer=screen_init,initargs=(str(r/'frozen/independent_screen.py'),)) if args.mode=='full' else contextlib.nullcontext(None)
    with poolctx as pool:
        futures=[pool.submit(screen_one,c) for c in cores] if pool else None
        for block in meta['blocks']:
            packed=(r/block['path']).read_bytes();assert len(packed)==block['bytes'] and digest(packed)==block['sha256'];number=0
            for line in gzip.decompress(packed).splitlines():
                z=json.loads(line);i=z['index'];assert i==len(seen);seen.append(i);number+=1
                pb=z['parent_utf8'].encode();qb=z['independent_utf8'].encode();assert digest(pb)==z['parent_sha256'] and digest(qb)==z['independent_sha256'];p,q=json.loads(pb),json.loads(qb)
                assert p['index']==q['index']==i and p['core']==q['core']==list(cores[i]) and p['status']==q['status']=='complete'
                assert p['binding']=={'source':parent_source,'prefix':meta['prefix_raw_sha256']} and q['binding']=={'source_sha256':ind_source,'prefix_sha256':meta['prefix_raw_sha256']}
                packedcode=base64.b64decode(z['codes_base64'],validate=True);assert digest(packedcode)==z['codes_gzip_sha256']==p['ledger_sha256'];codes=gzip.decompress(packedcode)
                assert len(codes)==105315 and digest(codes)==p['raw_code_sha256']==q['raw_code_sha256'];assert Counter(p['counts'])==counts_normalized(q['counts'],codes) and p['survivors']==q['survivors']
                if futures is not None:
                    cc,bb,ss=futures[i].result();assert cc==q['counts'] and bb==codes and ss==p['survivors']
                total.update(p['counts']);assert [v['row'] for v in z['records']]==p['survivors']
                for witness in z['records']:
                    row=tuple(witness['row']);assert row not in mother_seen;mother_seen.add(row);k=witness['death'];assert 1<=k<=67
                    tr=hgq.tree(witness['row'],depth=k,seconds=60);assert tr['status']=='dead_all_prefixes' and tr['depth_reached']==k and digest(canonical(tr['trace']))==witness['trace_sha256']
                    live=max(x['live'] for x in tr['trace']);assert live==witness['max_live'];maxlive=max(maxlive,live);maxdeath=max(maxdeath,k);layers+=k;depths[k]+=1
            assert number==block['count'];print(json.dumps({'verified_cores':len(seen),'verified_mothers':len(mother_seen),'mode':args.mode}),flush=True)
    assert seen==list(range(359)) and len(mother_seen)==9775 and layers==76172 and maxdeath==67 and maxlive==4
    assert total['raw']==37808085 and total['annulus_survived']==9775
    out={'status':'passed_complete_core_prefix_box_and_all_empty_trees','mode':args.mode,'screen_replayed':args.mode=='full','scope':'Full independently replays all 37,808,085 mother screens; trees binds archived screen bytes. Both fully regenerate 4913 cores,395 prefixes and all9775 empty trees.','core_box':4913,'prefixes':395,'cores':359,'mothers':len(mother_seen),'layers':layers,'max_death':maxdeath,'max_live':maxlive,'counts':dict(total),'death_histogram':dict(sorted(depths.items())),'workers':args.workers,'python':platform.python_version(),'interpreter':sys.executable,'argv':sys.argv,'elapsed_seconds':time.monotonic()-begin,'manifest_sha256':digest((r/'manifest.json').read_bytes()),'verifier_sha256':digest(Path(__file__).read_bytes())}
    args.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
