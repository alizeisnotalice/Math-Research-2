"""Freeze raw witness bytes and bind the independent multiplication replay."""
from argparse import ArgumentParser
import gzip
import hashlib
import json
from pathlib import Path


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def archive(raw,out):
    compressed=[]
    total=hashlib.sha256()
    with raw.open('rb') as src:
        k=0
        while content:=src.read(7_000_000):
            name=f'{raw.stem}-{k:02d}.part.gz'
            p=out/name
            with p.open('xb') as f:f.write(gzip.compress(content,mtime=0))
            assert p.stat().st_size<10*1024**2
            from_archive=gzip.decompress(p.read_bytes())
            assert content==from_archive
            total.update(from_archive)
            compressed.append({'file':name,'compressed_bytes':p.stat().st_size,'compressed_sha256':sha(p),
                               'raw_bytes':len(content),'raw_sha256':hashlib.sha256(content).hexdigest()})
            k+=1
    assert total.hexdigest()==sha(raw)
    return compressed


def main():
    if not __debug__:raise SystemExit('requires assertions')
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,default=Path.cwd())
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    census=[]
    sources={}
    sizes={29:(2347,18394),30:(5653,77977)}
    for j,(shapes,words) in sizes.items():
        raw=a.source/f'degree-{j}.json'
        check=a.source/f'degree-{j}-identity.json'
        original=json.loads(raw.read_bytes())
        replay=json.loads(check.read_bytes())
        assert original['mode']=='generate' and replay['mode']=='verify'
        assert original['offset']==replay['offset']==0 and original['limit'] is replay['limit'] is None
        assert not original['read_early_cache'] and not replay['read_early_cache']
        assert original['result']==replay['result'] and original['sources_sha256']==replay['sources_sha256']
        assert replay['witness_sha256']==sha(raw)
        for filename,h in original['sources_sha256'].items():
            assert sha(a.source/filename)==h
            assert sources.setdefault(filename,h)==h
        actual=original['result']
        details=actual['late_counts']
        assert details['shapes']==shapes and details['words']==words
        assert not actual['unresolved'] and not details.get('unresolved',0)
        certs=actual['certificates']
        assert len(certs)==shapes and all(not item['surviving_words'] for item in certs)
        projection=sum(len(p['projections']) for item in certs for p in item['fields'].values())
        units=sum(len(p['units']) for item in certs for p in item['fields'].values())
        assert projection==sum(details.get(f'mod{p}_projected_proofs',0) for p in (1279,163))
        assert units==sum(details.get(f'mod{p}_specialized_proofs',0) for p in (1279,163))
        census.append({'degree':j,'words':words,'projection_identities':projection,'unit_identities':units,
                       'original_sha256':sha(raw),'original_bytes':raw.stat().st_size,
                       'generation_seconds':original['seconds'],'identity_replay_seconds':replay['seconds'],
                       'identity_replay_sha256':sha(check),'chunks':archive(raw,a.out),
                       'early_counts':actual['early_counts'],'late_counts':details})
    report={'status':'verified','degrees':census,'sources_sha256':sources,
            'python':original['python'],'sympy':original['sympy'],
            'words':sum(r['words'] for r in census),
            'projection_identities':sum(r['projection_identities'] for r in census),
            'unit_identities':sum(r['unit_identities'] for r in census),'unresolved':0}
    with (a.out/'verification.json').open('x') as f:f.write(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:report[k] for k in ('status','words','projection_identities','unit_identities','unresolved')}))


if __name__=='__main__':main()
