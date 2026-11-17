#!/usr/bin/env python3
"""Aggregate content/coverage check; separate semantic replay remains required."""
import argparse,gzip,hashlib,itertools,json,math,os,time
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--ledger-root',type=Path,required=True);ap.add_argument('--record',type=Path,required=True);a=ap.parse_args()
assert os.path.ismount('/Volumes/data')
r=a.ledger_root;started=time.monotonic()
events=[json.loads(s) for s in (r/'events.jsonl').read_text().splitlines()]
ev=[e for e in events if e['event']=='chunk' and e['m']+e['n']==45]
seen=set();tot={k:0 for k in ['raw','canonical','excluded','fair','unknown']};pairs={};packedbytes=rawbytes=lines=0;finger=hashlib.sha256()
for e in sorted(ev,key=lambda e:(e['m'],e['index'])):
 m,n,idx=e['m'],e['n'],e['index'];key=(m,n,idx);assert key not in seen;seen.add(key)
 path=r/f'm{m:02d}n{n:02d}'/f'i{idx:06d}.jsonl.gz';packed=path.read_bytes();assert hashlib.sha256(packed).hexdigest()==e['gzip_sha256'];packedbytes+=len(packed)
 raw=gzip.decompress(packed);assert hashlib.sha256(raw).hexdigest()==e['raw_sha256'];rawbytes+=len(raw);k=raw.count(b'\n');assert k==e['counts']['canonical'];lines+=k
 assert e['counts']['raw']==sum(math.comb(n-1,j) for j in range(min(n-1,4)+1))
 for name in tot:tot[name]+=e['counts'][name]
 fingerprint=dict(m=m,n=n,index=idx,I=e['I'],counts=e['counts'],gzip_sha256=e['gzip_sha256'],raw_sha256=e['raw_sha256'])
 finger.update((json.dumps(fingerprint,sort_keys=True,separators=(',',':'))+'\n').encode());pairs.setdefault((m,n),[]).append(e)
assert set(pairs)=={(m,45-m) for m in range(1,23)}
for m in range(1,23):
 expected=[(0,*inside,m) for j in range(min(m-1,4)+1) for inside in itertools.combinations(range(1,m),j)]
 rows=pairs[m,45-m];assert len(rows)==len(expected)
 assert [(e['index'],tuple(e['I'])) for e in rows]==list(enumerate(expected))
assert tot==dict(raw=606304193,canonical=303168708,excluded=301344396,fair=1824312,unknown=0)
rec=dict(status='passed_parent_total45_full_decompression_hash_and_coverage',chunks=len(ev),counts=tot,packed_bytes=packedbytes,raw_bytes=rawbytes,lines=lines,binding_digest=finger.hexdigest(),elapsed_seconds=time.monotonic()-started,scope='Independent aggregate hash/coverage pass; full semantic per-J verification is supplied by audited_chunk_verifier and total45-report. No degree46 completion.',events_sha256=hashlib.sha256((r/'events.jsonl').read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
a.record.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec))
