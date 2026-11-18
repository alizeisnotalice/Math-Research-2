"""Preserve byte-exact sorted witnesses and all primary tree records in gzip parts."""
import argparse
from pathlib import Path
import hashlib
import gzip
import json


def split(source, destination, label):
    hashes=hashlib.sha256()
    parts=[]
    rows=0
    with source.open('rb') as stream:
        while True:
            raw=[];length=0
            for line in stream:
                raw.append(line)
                length+=len(line)
                if length>=7_000_000:break
            if not raw:break
            name=f'{label}-{len(parts):02}.jsonl.gz'
            payload=b''.join(raw)
            hashes.update(payload)
            with (destination/name).open('xb') as f, gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as z:
                z.write(payload)
            assert (destination/name).stat().st_size<10*1024**2
            parts.append({'file':name,'rows':len(raw),'raw_sha256':hashlib.sha256(payload).hexdigest()})
            rows+=len(raw)
    return {'rows':rows,'parts':parts,'raw_sha256':hashes.hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--roots',type=Path,required=True)
    parser.add_argument('--trees',type=Path,required=True)
    parser.add_argument('--archive',type=Path,required=True)
    args=parser.parse_args()
    args.archive.mkdir()
    index={key:split(path,args.archive,key) for key,path in [('roots',args.roots),('trees',args.trees)]}
    assert index['roots']['rows']==index['trees']['rows']==869396
    with (args.archive/'index.json').open('x') as f:
        f.write(json.dumps(index,sort_keys=True,indent=2)+'\n')


if __name__=='__main__':
    main()
