"""Bind completed exhaustive verification receipts to current ledger bytes and coverage.

This is an integrity/coverage aggregation, not a replacement for the independent
per-record exhaustive verifier that produced verification.json for each shard.
"""
import argparse, gzip, hashlib, json, math, time
from collections import Counter
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--ledger-root',type=Path,required=True)
    p.add_argument('--record',type=Path,required=True);a=p.parse_args();started=time.monotonic()
    source_hashes={name:sha(SOURCE/name) for name in ['support_search_layer.py','verify_support_layer.py']}
    expected={(m,total-m) for total in [41,42] for m in range(1,total//2+1)}
    seen=set();totals=Counter();rows=[]
    for m in range(1,22):
        name=f'six42-m{m:02d}';folder=HERE/name;live=a.ledger_root/name
        s=json.loads((folder/'summary.json').read_text());v=json.loads((folder/'verification.json').read_text())
        j=json.loads((folder/'job.json').read_text());c=s['config']
        assert c['total']==42 and c['min_total']==41 and c['max_terms']==6
        assert c['min_degree']==c['max_first_degree']==m and c['version']==3
        assert s['complete'] and v['verified'] and j['status']=='verified' and j['process_exited']
        assert s['counts']==v['counts'] and s['counts']['unknown']==0
        assert (folder/'unknown.jsonl').read_bytes()==b''
        assert s['script_sha256']==source_hashes['support_search_layer.py']
        assert v['verifier_sha256']==source_hashes['verify_support_layer.py']
        for shard in s['shards']:
            pair=(shard['m'],shard['n']);assert pair in expected and pair not in seen and shard['complete']
            seen.add(pair)
            number=lambda d:sum(math.comb(d-1,k) for k in range(min(4,d-1)+1))
            assert shard['raw']==number(pair[0])*number(pair[1])
        raw_hash=hashlib.sha256();raw_bytes=0
        with gzip.open(live/'ledger.jsonl.gz','rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''):
                raw_hash.update(block);raw_bytes+=len(block)
        assert raw_hash.hexdigest()==s['ledger_sha256_uncompressed']==v['ledger_sha256_uncompressed']
        totals.update(s['counts'])
        rows.append({'shard':name,'ledger_path':str(live/'ledger.jsonl.gz'),
                     'gzip_sha256':sha(live/'ledger.jsonl.gz'),'uncompressed_sha256':raw_hash.hexdigest(),
                     'uncompressed_bytes':raw_bytes,'counts':s['counts']})
    assert seen==expected and totals['raw']==585833432 and totals['canonical']==283126617
    assert totals['unknown']==0 and totals['canonical']==totals['fair']+totals['excluded']
    out={'status':'completed_exhaustive_receipts_bound_to_full_ledger_bytes',
         'degree_totals':[41,42],'max_terms_each':6,'degree_pairs':len(seen),'shards':len(rows),
         'source_hashes':source_hashes,'counts':dict(totals),'ledgers':rows,
         'elapsed_seconds':time.monotonic()-started,
         'scope':'No new per-case reclassification here; independent exhaustive verifier completed every shard. This additionally checks complete degree coverage and unchanged ledger bytes.'}
    a.record.write_text(json.dumps(out,indent=2)+'\n');print({k:out[k] for k in ['status','counts','shards','elapsed_seconds']})
if __name__=='__main__':main()
