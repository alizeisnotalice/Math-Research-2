"""Parent audit: lossless archive bindings and independent support coverage."""
import argparse,gzip,hashlib,itertools,json,sys,time
from collections import Counter
from pathlib import Path

def main():
    if not __debug__:raise RuntimeError('no -O')
    p=argparse.ArgumentParser();p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.package;t=time.monotonic();sha=lambda b:hashlib.sha256(b).hexdigest()
    assert not a.out.exists(), 'fresh output required'
    m=json.loads((r/'package-manifest.json').read_text())
    shared=(r/m['shared_package']).resolve()
    shared_names={x['source_name'] for x in m['shared_input_archives']}
    def read_input(name):
        root=shared if name in shared_names else r
        return gzip.decompress((root/'input'/f'{name}.gz').read_bytes())
    for x in m['shared_input_archives']:
        b=(shared/x['file']).read_bytes();assert len(b)==x['bytes'] and sha(b)==x['sha256'] and sha(gzip.decompress(b))==x['source_sha256']
    for x in m['files']:
        b=(r/x['path']).read_bytes();assert len(b)==x['bytes'] and sha(b)==x['sha256']
    for x in m['input_archives']:
        assert sha(gzip.decompress((r/x['file']).read_bytes()))==x['source_sha256']
    expected={}
    for degree in range(1,27):
        i=0
        for size in range(2,min(degree+1,6)+1):
            for middle in itertools.combinations(range(1,degree),size-2):
                expected[(degree,i)]=(0,*middle,degree);i+=1
    assert len(expected)==83681
    events={}
    for line in read_input('run__events.jsonl').decode().splitlines():
        e=json.loads(line)
        if e.get('event')!='chunk' or e.get('total')!=53:continue
        key=(e['m'],e['index']);assert key in expected and key not in events
        assert e['n']==53-e['m'] and tuple(e['I'])==expected[key]
        events[key]=e
    assert set(events)==set(expected)
    progress=json.loads(read_input('progress53.json'));binding=progress['binding']
    assert progress['status']=='complete' and progress['verified']==83681 and progress['remaining']==0 and progress['failures']==[]
    names={'events':'run__events.jsonl','launch':'run__launch-manifest.json','payload':'run__payload.py','plan':'plan__support54-plan.json','runner':'code__runner.py','source':'plan__frozen_support.py','terminal':'run__exit-receipt.json','verifier':'code__audited_chunk_verifier.py'}
    assert binding['total']==53 and binding['cap']==6
    assert binding['runner']=='7d19ae57220d2d06b0f34c251cbb3decdeae286c3d44ea84bab0da06027a63eb'
    assert binding['verifier']=='a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b'
    assert binding['source']=='d68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710'
    for k,n in names.items():assert binding[k]==sha(read_input(n))
    plan=json.loads(read_input('plan__support54-plan.json'))
    actual={}
    for j in plan['jobs']:
        if j['total']!=53:continue
        key=j['m'],j['index'];assert key not in actual and j['n']==53-j['m'];actual[key]=tuple(j['I'])
    assert actual==expected and len(plan['jobs'])==185264
    original=json.loads(read_input('run__exit-receipt.json'))
    assert original['returncode']==255 and not original['selected_set_complete'] and original['chunks_persisted']==94465
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
    assert seen==set(expected) and rows==83681
    assert counts['unknown']==0
    terminal=json.loads(read_input('exit53.json'))
    assert terminal['returncode']==0 and terminal['child_waited'] is True and terminal['verified']==83681 and terminal['remaining']==0
    launch53=json.loads(read_input('launch53.json'))
    assert terminal['status']=='complete' and terminal['pid']==launch53['pid']==progress['pid']
    assert launch53['runner_sha256']==binding['runner']
    raw_audit=json.loads((r.parent/'audits/raw53.json').read_text())
    assert raw_audit['status']=='passed' and raw_audit['chunks']==83681
    assert raw_audit['binding']['events']==binding['events']
    assert raw_audit['counts']==dict(counts)
    assert raw_audit['gzip_bytes']==5420857949 and raw_audit['raw_bytes']==78491726768
    assert sum(e['counts']['canonical'] for e in events.values())==counts['canonical']
    out={'status':'passed_parent_full_archive_hash_binding_and_support_coverage','scope':'83681 independently generated I shards; semantic replay supplied by byte-bound remote verifier, original bytes additionally read by raw53 audit.','records':rows,'aggregate':dict(counts),'raw_audit_sha256':sha((r.parent/'audits/raw53.json').read_bytes()),'manifest_sha256':sha((r/'package-manifest.json').read_bytes()),'verifier_sha256':sha(Path(__file__).read_bytes()),'elapsed_seconds':time.monotonic()-t,'python':sys.version}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    tmp=a.out.with_suffix('.partial');tmp.write_text(json.dumps(out,indent=2)+'\n');tmp.replace(a.out);print(json.dumps(out))

if __name__=='__main__':main()
