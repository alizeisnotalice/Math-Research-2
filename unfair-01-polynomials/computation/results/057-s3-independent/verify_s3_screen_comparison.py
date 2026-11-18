"""Check all returned integer-enclosure passes against the production stream."""
import argparse
import hashlib
import gzip
import heapq
import json
from pathlib import Path
from contextlib import ExitStack


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(file, reorder=False):
    previous=None
    rows=(tuple(json.loads(line)) for line in file)
    for row in sorted(rows) if reorder else rows:
        assert len(row)==9 and (previous is None or row>previous)
        yield row
        previous=row


def records_unordered(raw):
    for line in raw.splitlines():
        row=tuple(json.loads(line))
        assert len(row)==9
        yield row


def main():
    if not __debug__:
        raise RuntimeError('run without -O')
    p=argparse.ArgumentParser()
    p.add_argument('--production',type=Path,required=True)
    p.add_argument('--return-dir',type=Path,required=True)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True)
    a=p.parse_args()
    result=json.loads((a.return_dir/'summary.json').read_text())
    parts=json.loads(a.manifest.read_text())['parts']
    assert result['status']=='complete' and result['child_processes_waited'] and len(parts)==13
    assert result['manifest_sha256']==sha(a.manifest) and result['source_sha256']==sha(a.source)
    assert [part['input'] for part in result['parts']]==parts
    total=result['counts']
    assert total['unresolved']==0 and total['pass']==869396 and total['mothers']==1523605945
    assert sha(a.production)=='6da30bcf5f93f0fc121923d1979a37d0508327a5401d10c4e8d86ceb41fd546a'
    with ExitStack() as stack, a.production.open('rb') as production:
        streams=[]
        for part in result['parts']:
            f=a.return_dir/part['output']
            if f.exists():
                raw=f.read_bytes()
            else:
                raw=gzip.decompress(f.with_suffix(f.suffix+'.gz').read_bytes())
            assert hashlib.sha256(raw).hexdigest()==part['output_sha256'] and part['counts']['unresolved']==0
            streams.append(iter(sorted(records_unordered(raw))))
        count=0
        for x,y in zip(records(production),heapq.merge(*streams),strict=True):
            assert x==y
            count+=1
        assert count==total['pass']
    print(json.dumps({'status':'all_integer_enclosure_and_production_pass_rows_equal', 'rows':count,'mothers':total['mothers']}))


if __name__=='__main__':
    main()
