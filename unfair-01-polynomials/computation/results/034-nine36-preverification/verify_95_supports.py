"""Exact support certificates, separately bound to unverified-search snapshots."""
import argparse,hashlib,importlib.util,json,subprocess,sys
from pathlib import Path
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(r):return Path(r['source_unknown']['path']).parent.name,r['source_unknown']['line']
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--record',type=Path,required=True);a=ap.parse_args()
    data=json.loads((a.bundle/'nine36-direct-preverification.json').read_text());rows=data['rows']
    original=json.loads((a.bundle/'nine36-preverification-map.json').read_text());assert data['input_sha256']==sha(a.bundle/'nine36-preverification-map.json')
    assert len(rows)==len({key(r) for r in rows})==len({(tuple(r['I']),tuple(r['J'])) for r in rows})==95
    assert {key(r):(r['I'],r['J']) for r in rows}=={key(r):(r['I'],r['J']) for r in original['rows']}
    expected={};sources={}
    for shard in ['nine36-m15','nine36-m16','nine36-m17']:
        d=a.bundle/'inputs'/shard;s=json.loads((d/'summary.json').read_text());assert s['complete']
        sources[shard]={name:sha(d/name) for name in ['summary.json','unknown.jsonl','job.json']}
        sources[shard]['declared_ledger_sha256_uncompressed']=s['ledger_sha256_uncompressed']
        unknown=[json.loads(line) for line in (d/'unknown.jsonl').read_text().splitlines()]
        assert len(unknown)==s['counts']['unknown']
        for line,row in enumerate(unknown,1):expected[shard,line]=(row['I'],row['J'])
    assert expected=={key(r):(r['I'],r['J']) for r in rows}
    ddpath=a.topic/'computation/diagonal_domination.py';assert sha(ddpath)==data['source_sha256']
    spec=importlib.util.spec_from_file_location('saved_diagonal_verifier',ddpath);dd=importlib.util.module_from_spec(spec);spec.loader.exec_module(dd)
    positive={('nine36-m16',8):a.bundle/'line8-recovered-candidate.json',('nine36-m16',15):a.bundle/'positive-preverification/nine36-m16-line15.json',('nine36-m16',16):a.bundle/'positive-preverification/nine36-m16-line16.json'}
    verifier=a.topic/'computation/verify_positive_certificate.py';direct=0;receipts=[]
    for r in rows:
        assert r['search_summary_sha256']==sources[key(r)[0]]['declared_ledger_sha256_uncompressed']
        c=r['direct_diagonal_certificate']
        if c:
            assert dd.verify_certificate(r['I'],r['J'],c);direct+=1
        else:
            path=positive.pop(key(r));c=json.loads(path.read_text());assert c['case']['I']==r['I'] and c['case']['J']==r['J']
            run=subprocess.run([sys.executable,str(verifier),str(path)],text=True,capture_output=True,check=True);v=json.loads(run.stdout);assert v['status']=='verified'
            receipts.append({'support_row':list(key(r)),'certificate':str(path.relative_to(a.bundle)),'sha256':sha(path),'verification':v})
    assert direct==92 and len(receipts)==3 and not positive
    out={'status':'95_support_local_certificates_verified','direct':direct,'positive':len(receipts),'layer_completion_claim':False,'source_hashes':sources,'historical_field_semantics':{'search_summary_sha256':'This discovery-field name actually contains summary.ledger_sha256_uncompressed, not the summary-file SHA. Both are explicitly bound here.'},'diagonal_verifier_sha256':sha(ddpath),'positive_verifier_sha256':sha(verifier),'receipts':receipts,'scope':'The three complete search snapshots were awaiting independent full-ledger verification. This certifies these supports only, not complete shards or the whole 44/36 layer.'}
    a.record.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','direct','positive','layer_completion_claim']}))
if __name__=='__main__':main()
