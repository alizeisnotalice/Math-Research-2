"""Partition the published group stream into disjoint, ordered subsequences."""
import argparse
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path

EXPECTED = '5cf9a01d716e19c9bf5dac8bb43e90ce49c4d238a09b434ce89d09d85528ca36'


def main():
    p=argparse.ArgumentParser()
    p.add_argument('groups',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    if not __debug__:
        raise RuntimeError('run without -O')
    args.output.mkdir()
    digest=hashlib.sha256()
    size=0
    names=[args.output/f'groups-{i:02}.txt' for i in range(13)]
    with ExitStack() as stack, args.groups.open('rb') as source:
        shards=[stack.enter_context(name.open('xb')) for name in names]
        for size,line in enumerate(source,1):
            digest.update(line)
            shards[(size-1)%13].write(line)
    assert size==1915327 and digest.hexdigest()==EXPECTED
    parts=[{'file':f.name,'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size} for f in names]
    manifest={'input_sha256':EXPECTED,'groups':size,'parts':parts}
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    main()
