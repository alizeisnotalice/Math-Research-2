"""Stream the full AB ledger, compare every HGQ layer, emit small certificates."""
import argparse,hashlib,itertools,json,re,time,os
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def stream_cores(path):
    dec=json.JSONDecoder()
    with path.open() as f:
        buf=''
        while True:
            buf+=f.read(65536);m=re.search(r'"cores"\s*:\s*\[',buf)
            if m:break
            if len(buf)>1<<20:raise ValueError('missing cores array')
        header=json.loads(buf[:m.end()]+']}')
        assert header['distinct_mothers']==90092 and header['duplicate_trace_conflicts']==0
        assert header['parent_current_uncovered_mothers']==0 and header['parent_hash_mismatches']==[]
        buf=buf[m.end():];first=True
        while True:
            buf=buf.lstrip()
            if not buf:
                more=f.read(65536)
                if not more:raise ValueError('truncated array')
                buf+=more;continue
            if buf.startswith(']'):
                assert (buf[1:]+f.read()).strip()=='}';return
            if not first:
                assert buf[0]==',';buf=buf[1:].lstrip()
            while True:
                buf=buf.lstrip()
                try:obj,end=dec.raw_decode(buf);break
                except json.JSONDecodeError:
                    if len(buf)>1<<20:raise ValueError('malformed or oversized core record')
                    more=f.read(65536)
                    if not more:raise
                    buf+=more
            yield obj;buf=buf[end:];first=False

def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser()
    for n in ['ab-ledger','hgq','screen','extension','out']:ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);start=time.monotonic();ab_hash=sha(a.ab_ledger);extension=json.loads(a.extension.read_text());ab_iter=iter(stream_cores(a.ab_ledger))
    total=0;max_depth=0;max_live=0;layers=0;tn=[];last_row=None;full_digest=hashlib.sha256()
    for i in range(669):
        src=a.screen/f'tn-{i:03}.json';hpath=a.hgq/f'tn-{i:03}.json';s=json.loads(src.read_text());h=json.loads(hpath.read_text());src_hash=sha(src)
        assert s['status']=='complete' and h['index']==s['index']==i and h['screen_sha256']==src_hash
        assert [r['row'] for r in h['records']]==s['survivors']
        rows=[]
        for rec in h['records']:
            ab=next(ab_iter);key=tuple(rec['row'])
            assert last_row is None or key>last_row;last_row=key
            assert ab['mother']==rec['row'] and ab['source']['tn_index']==i and ab['source']['sha256']==src_hash
            if rec['status']=='unknown_depth':
                assert rec['row']==extension['row'] and rec['trace']==extension['trace'][:len(rec['trace'])];rec=extension
            assert ab['status']==rec['status']=='dead_all_prefixes'
            assert ab['depth_reached']==rec['depth_reached'] and ab['trace']==rec['trace']
            trace=rec['trace'];assert [z['k'] for z in trace]==list(range(1,len(trace)+1)) and trace[-1]['live']==0
            mh=max(z['live'] for z in trace);digest=hashlib.sha256(canonical(trace)).hexdigest()
            record={'row':list(key),'death':len(trace),'max_live':mh,'trace_sha256':digest};rows.append(record);full_digest.update(canonical(record)+b'\n')
            max_depth=max(max_depth,len(trace));max_live=max(max_live,mh);layers+=len(trace);total+=1
        out={'index':i,'screen_sha256':src_hash,'hgq_sha256':sha(hpath),'ab_ledger_sha256':ab_hash,'mothers':len(rows),'records':rows}
        target=a.out/f'tn-{i:03}.json';target.write_text(json.dumps(out,separators=(',',':'))+'\n')
        tn.append({'index':i,'mothers':len(rows),'certificate_sha256':sha(target),'screen_sha256':src_hash,'hgq_sha256':sha(hpath)})
    assert next(ab_iter,None) is None and total==90092
    out={'status':'passed_all_90092_mothers_every_AB_HGQ_layer_and_source_binding','mothers':total,'TN':669,'trace_layers':layers,'max_death':max_depth,'max_live':max_live,'ordered_compact_certificates_sha256':full_digest.hexdigest(),'ab_ledger_sha256':ab_hash,'extension_sha256':sha(a.extension),'per_TN':tn,'elapsed_seconds':time.monotonic()-start,'script_sha256':sha(Path(__file__))}
    (a.out/'comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','mothers','trace_layers','max_death','max_live','elapsed_seconds']}))

if __name__=='__main__':main()
