"""Parent audit: lossless archive bindings and independent support coverage."""
import argparse,gzip,hashlib,itertools,json,sys,time
from collections import Counter
from pathlib import Path

def main():
    if not __debug__:raise RuntimeError('no -O')
    p=argparse.ArgumentParser();p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.package;t=time.monotonic();sha=lambda b:hashlib.sha256(b).hexdigest()
    m=json.loads((r/'package-manifest.json').read_text())
    for x in m['files']:
        b=(r/x['path']).read_bytes();assert len(b)==x['bytes'] and sha(b)==x['sha256']
    for x in m['input_archives']:
        assert sha(gzip.decompress((r/x['file']).read_bytes()))==x['source_sha256']
    expected={}
    for degree in range(1,26):
        i=0
        for size in range(2,min(degree+1,6)+1):
            for middle in itertools.combinations(range(1,degree),size-2):
                expected[(degree,i)]=(0,*middle,degree);i+=1
    assert len(expected)==68405
    events={}
    with gzip.open(r/'input/run__events.jsonl.gz','rt') as f:
        for line in f:
            e=json.loads(line)
            if e.get('event')!='chunk' or e.get('total')!=51:continue
            key=(e['m'],e['index']);assert key in expected and key not in events
            assert e['n']==51-e['m'] and tuple(e['I'])==expected[key]
            events[key]=e
    assert set(events)==set(expected)
    progress=json.loads(gzip.decompress((r/'input/progress51.json.gz').read_bytes()));binding=progress['binding']
    assert progress['status']=='complete' and progress['verified']==68405 and progress['remaining']==0 and progress['failures']==[]
    names={'events':'run__events.jsonl','launch':'run__launch-manifest.json','payload':'run__payload.py','plan':'plan__support52-plan.json','runner':'code__runner.py','source':'plan__frozen_support.py','terminal':'run__exit-receipt.json','verifier':'code__audited_chunk_verifier.py'}
    assert binding['total']==51 and binding['cap']==6
    for k,n in names.items():assert binding[k]==sha(gzip.decompress((r/'input'/f'{n}.gz').read_bytes()))
    seen=set();counts=Counter();rows=0
    for archive in m['coverage']['cache_archives']:
        number=0
        with gzip.open(r/archive['file'],'rt') as f:
            for line in f:
                z=json.loads(line);raw=z['utf8'].encode();assert sha(raw)==z['sha256']
                rel=Path(z['relativepath']);assert not rel.is_absolute() and '..' not in rel.parts
                key=(int(rel.parent.name[1:]),int(rel.stem[1:]));assert key in expected and key not in seen
                seen.add(key);q=json.loads(raw);e=events[key]
                assert q['binding']==binding and q['status']=='passed' and q['errors']==[]
                assert q['event_sha256']==sha(json.dumps(e,sort_keys=True).encode())
                result=q['result'];assert result['counts']==e['counts'] and result['verified_rows']==e['counts']['canonical'] and result['raw_support_count']==e['counts']['raw']
                counts.update(result['counts']);number+=1;rows+=1
        assert number==archive['records']
    assert seen==set(expected) and rows==68405
    assert counts['unknown']==0
    terminal=json.loads(gzip.decompress((r/'input/exit51.json.gz').read_bytes()))
    assert terminal['returncode']==0 and terminal['child_waited'] is True and terminal['verified']==68405 and terminal['remaining']==0
    out={'status':'passed_parent_full_archive_hash_binding_and_support_coverage','scope':'68405 independently generated I shards; semantic replay supplied by byte-bound remote verifier, original bytes additionally read by raw51 audit.','records':rows,'aggregate':dict(counts),'manifest_sha256':sha((r/'package-manifest.json').read_bytes()),'verifier_sha256':sha(Path(__file__).read_bytes()),'elapsed_seconds':time.monotonic()-t,'python':sys.version}
    a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
