#!/usr/bin/env python3
"""Independent aggregate pass over all old+new total46 ledger bytes and I coverage."""
import argparse, gzip, hashlib, itertools, json, math, time
from pathlib import Path
ap=argparse.ArgumentParser()
for name in ['old-run','new-run','semantic-report','out']:ap.add_argument('--'+name,type=Path,required=True)
a=ap.parse_args();start=time.monotonic();entries=[];eventhash={}
for source,root in [('old',a.old_run),('new',a.new_run)]:
 data=(root/'events.jsonl').read_bytes();eventhash[source]=hashlib.sha256(data).hexdigest()
 for line in data.splitlines():
  e=json.loads(line)
  if e.get('event')=='chunk' and e['m']+e['n']==46:entries.append((e,source,root))
report=json.loads(a.semantic_report.read_text());assert report['status']=='passed_complete_combined_total46' and not report['manifest']['errors']
semantic={(x['source'],x['job']):x for x in report['chunk_results']};assert len(semantic)==len(report['chunk_results'])==44551
seen=set();by_pair={};counts=dict(raw=0,canonical=0,excluded=0,fair=0,unknown=0);packedbytes=rawbytes=lines=0;finger=hashlib.sha256()
for e,source,root in sorted(entries,key=lambda t:(t[0]['m'],t[0]['index'])):
 m,n,idx=e['m'],e['n'],e['index'];key=(m,n,idx);assert key not in seen;seen.add(key)
 job=f'm{m:02d}n{n:02d}/i{idx:06d}';v=semantic[(source,job)];assert not v['errors'] and v['counts']==e['counts']
 assert v['event_gzip_sha256']==e['gzip_sha256'] and v['event_raw_sha256']==e['raw_sha256']
 packed=(root/(job+'.jsonl.gz')).read_bytes();assert hashlib.sha256(packed).hexdigest()==e['gzip_sha256'];raw=gzip.decompress(packed);assert hashlib.sha256(raw).hexdigest()==e['raw_sha256']
 k=raw.count(b'\n');assert k==e['counts']['canonical']==v['verified_rows'];packedbytes+=len(packed);rawbytes+=len(raw);lines+=k
 assert e['counts']['raw']==sum(math.comb(n-1,j) for j in range(min(4,n-1)+1))
 for field in counts:counts[field]+=e['counts'][field]
 by_pair.setdefault((m,n),[]).append((idx,tuple(e['I'])))
 finger.update((json.dumps(dict(key=key,source=source,I=e['I'],counts=e['counts'],gzip=e['gzip_sha256'],raw=e['raw_sha256']),sort_keys=True,separators=(',',':'))+'\n').encode())
assert set(by_pair)=={(m,46-m) for m in range(1,24)}
for m in range(1,24):
 expected=[(0,*inside,m) for j in range(min(m-1,4)+1) for inside in itertools.combinations(range(1,m),j)]
 assert by_pair[m,46-m]==list(enumerate(expected))
assert counts==report['counts']==dict(raw=787855938,canonical=372691286,excluded=370540360,fair=2150926,unknown=0)
assert report['cache']['binding']['old_events_sha256']==eventhash['old'] and report['cache']['binding']['new_events_sha256']==eventhash['new']
out=dict(status='passed_parent_full_combined46_hash_coverage_semantic_binding',chunks=len(seen),counts=counts,packed_bytes=packedbytes,raw_bytes=rawbytes,lines=lines,binding_digest=finger.hexdigest(),events_sha256=eventhash,semantic_report_sha256=hashlib.sha256(a.semantic_report.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-start,scope='All total46 cap6 I chunks; independent full-byte aggregate check bound to a separate complete semantic replay. No total47/48 or general conclusion.')
a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
